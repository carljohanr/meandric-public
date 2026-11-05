#include <cmath>
#include <cstdint>
#include <vector>
#include <algorithm>
#include <limits>

// Streaming paired forward-row evaluator. Search only: certificates are checked
// by an independent exact implementation. No row matrix or thread pool.
namespace {
uint64_t reverse(uint64_t a, int len) {
    a=((a>>1)&0x5555555555555555ULL)|((a&0x5555555555555555ULL)<<1);
    a=((a>>2)&0x3333333333333333ULL)|((a&0x3333333333333333ULL)<<2);
    a=((a>>4)&0x0f0f0f0f0f0f0f0fULL)|((a&0x0f0f0f0f0f0f0f0fULL)<<4);
    a=((a>>8)&0x00ff00ff00ff00ffULL)|((a&0x00ff00ff00ff00ffULL)<<8);
    a=((a>>16)&0x0000ffff0000ffffULL)|((a&0x0000ffff0000ffffULL)<<16);
    a=(a>>32)|(a<<32);
    return a>>(64-len);
}
struct Feature { int a[40]; int size=0; double value=0; };
struct Model {
    int r,k,q,nu; uint64_t mask, fmask; const int32_t* lookup; const int32_t* edges; const double* x;
    Feature feature(uint64_t s,int L,int h) const {
        Feature f;
        f.a[f.size++]=lookup[(s>>(L-h-r-1))&mask];
        uint64_t rev=reverse(s,L);
        for(int j=0;j<=h-k;j++) f.a[f.size++]=nu+((s>>(L-j-k))&fmask);
        for(int j=0;j<=L-h-1-k;j++) f.a[f.size++]=nu+((rev>>(L-j-k))&fmask);
        for(int j=0;j<f.size;j++) f.value+=x[f.a[j]];
        return f;
    }
    double forward(uint64_t s,uint64_t physical,double multiplier,double* gradient) const {
        int L=2*r+3,h=r+1,b=r+1;
        Feature base=feature(s,L,h);
        double sum=0;
        auto move=[&](uint64_t t,int len,int head,int correction=0) {
            Feature dest=feature(t,len,head);
            int variable=nu+(1<<k)+std::abs(correction)-1;
            double eta=correction==0?0:(correction>0?x[variable]:-x[variable]);
            double term=std::exp(dest.value-base.value+eta);sum+=term;
            if(gradient && multiplier!=0) {
                double p=multiplier*term;
                for(int j=0;j<dest.size;j++) gradient[dest.a[j]]+=p;
                for(int j=0;j<base.size;j++) gradient[base.a[j]]-=p;
                if(correction) gradient[variable]+=(correction>0?p:-p);
            }
        };
        uint64_t cmask=(1ULL<<(2*q))-1;
        int left=-edges[(physical>>(L-h-q))&cmask];
        int right=edges[(physical>>(L-h-q-1))&cmask];
        move(s,L,h-1,left);move(s,L,h+1,right);
        uint64_t upper=s>>(b+1), lower=s&((1ULL<<b)-1);
        if((s>>b)&1) move((upper<<(b+3))|(5ULL<<b)|lower,L+2,h+1);
        else {
            move((upper<<(b+3))|(4ULL<<b)|lower,L+2,h+1);
            move((upper<<(b+3))|(1ULL<<b)|lower,L+2,h+1);
        }
        int tri=(s>>(b-1))&7;
        if(tri==2 || tri==3 || tri==6) {
            uint64_t t=((s>>(b+2))<<b)|((uint64_t)(tri!=2)<<(b-1))|(s&((1ULL<<(b-1))-1));
            move(t,L-2,h-1);
        }
        return sum;
    }
};
}
extern "C" double evaluate(int r,int k,int q,int nu,int ne,const int32_t* lookup,const int32_t* edges,const double* x,
                           double beta,double* gradient,double* maximum,double* rowsout) {
    Model model{r,k,q,nu,(1ULL<<(2*r+1))-1,(1ULL<<k)-1,lookup,edges,x};
    uint64_t groups=1ULL<<(2*r+2),complement=2*groups-1;
    std::vector<double> logs(groups);
    double biggest=-std::numeric_limits<double>::infinity();
    for(uint64_t z=0;z<groups;z++) {
        // Both weight components use the ORIGINAL physical head-edge field.
        double row=(model.forward(z,z,0,nullptr)+model.forward(complement-z,z,0,nullptr))/2;
        logs[z]=std::log(row);biggest=std::max(biggest,logs[z]);
        if(rowsout) rowsout[z]=row;
    }
    *maximum=std::exp(biggest);
    if(!gradient) return biggest;
    std::fill(gradient,gradient+nu+(1<<k)+ne,0.0);
    double total=0;
    for(uint64_t z=0;z<groups;z++) {
        double probability=std::exp(beta*(logs[z]-biggest));
        total+=probability;
        if(probability<1e-25) continue;
        double p=probability*std::exp(-logs[z])/2;
        model.forward(z,z,p,gradient);model.forward(complement-z,z,p,gradient);
    }
    for(int j=0;j<nu+(1<<k)+ne;j++) gradient[j]/=total;
    return biggest+std::log(total)/beta;
}
