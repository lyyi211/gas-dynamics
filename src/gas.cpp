#include <algorithm>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>
#include <omp.h>

struct V { double x=0,y=0; };
struct Gas {
    int n, cells;
    double L, rc=2.5, uc, fc;
    std::vector<V> r,v,a,unwrapped;
    double potential=0, virial=0;
    Gas(int count,double density,double temperature,unsigned seed):n(count),L(std::sqrt(n/density)),r(n),v(n),a(n),unwrapped(n) {
        if(L<=2*rc) throw std::runtime_error("Box must be wider than twice the cutoff; increase N or reduce density");
        cells=std::max(1,int(L/rc));
        double q=std::pow(rc,-6);
        uc=4*(q*q-q); fc=24*(2*q*q-q)/rc;
        std::mt19937 gen(seed); std::normal_distribution<double> normal;
        int side=int(std::ceil(std::sqrt(n))); V mean;
        for(int i=0;i<n;i++) {
            r[i]={(i%side+0.5)*L/side,(i/side+0.5)*L/side};
            unwrapped[i]=r[i];
            v[i]={normal(gen),normal(gen)}; mean.x+=v[i].x/n; mean.y+=v[i].y/n;
        }
        double sum=0;
        for(auto &w:v){w.x-=mean.x;w.y-=mean.y;sum+=w.x*w.x+w.y*w.y;}
        double scale=std::sqrt((2*n-2)*temperature/sum);
        for(auto &w:v){w.x*=scale;w.y*=scale;}
    }
    double kinetic() const {double k=0;for(auto w:v)k+=0.5*(w.x*w.x+w.y*w.y);return k;}
    void forces(const std::string &method) {
        std::vector<std::vector<int>> bins(cells*cells);
        if(method=="grid")for(int i=0;i<n;i++)bins[int(r[i].y/L*cells)*cells+int(r[i].x/L*cells)].push_back(i);
        double energy=0, pair_virial=0;
        #pragma omp parallel for reduction(+:energy,pair_virial) if(method=="openmp") schedule(static)
        for(int i=0;i<n;i++) {
            V force; double e=0, w=0;
            auto pair=[&](int j){
                if(i==j)return;
                double dx=r[i].x-r[j].x,dy=r[i].y-r[j].y;
                dx-=L*std::round(dx/L);dy-=L*std::round(dy/L);
                double d2=dx*dx+dy*dy;
                if(d2>=rc*rc)return;
                double d=std::sqrt(d2),q=1/(d2*d2*d2);
                double f=(24*(2*q*q-q)/d-fc)/d;
                force.x+=f*dx;force.y+=f*dy;w+=0.5*f*d2;
                e+=0.5*(4*(q*q-q)-uc+(d-rc)*fc);
            };
            if(method=="grid") {
                int cx=int(r[i].x/L*cells),cy=int(r[i].y/L*cells);
                int visited[9],used=0;
                for(int oy=-1;oy<=1;oy++)for(int ox=-1;ox<=1;ox++) {
                    int b=((cy+oy+cells)%cells)*cells+(cx+ox+cells)%cells;
                    if(std::find(visited,visited+used,b)!=visited+used)continue;
                    visited[used++]=b;for(int j:bins[b])pair(j);
                }
            } else for(int j=0;j<n;j++)pair(j);
            a[i]=force;energy+=e;pair_virial+=w;
        }
        potential=energy;virial=pair_virial;
    }
    void step(double dt,const std::string &method) {
        for(int i=0;i<n;i++) {
            v[i].x+=0.5*dt*a[i].x;v[i].y+=0.5*dt*a[i].y;
            r[i].x+=dt*v[i].x;r[i].y+=dt*v[i].y;
            unwrapped[i].x+=dt*v[i].x;unwrapped[i].y+=dt*v[i].y;
            r[i].x-=L*std::floor(r[i].x/L);r[i].y-=L*std::floor(r[i].y/L);
        }
        forces(method);
        for(int i=0;i<n;i++){v[i].x+=0.5*dt*a[i].x;v[i].y+=0.5*dt*a[i].y;}
    }
};
int main(int argc,char **argv) {try {
    int n=256,steps=3000,stride=15,threads=4; unsigned seed=42;
    double density=.3,temp=1.,dt=.002;std::string method="grid",out="output";
    for(int i=1;i<argc;i++) {
        std::string key=argv[i];
        if(key=="--help"){std::cout<<"--n 256 --steps 3000 --stride 15 --density .3 --temperature 1 --dt .002 --method sequential|openmp|grid --threads 4 --seed 42 --out output\n";return 0;}
        if(i+1==argc)throw std::runtime_error("Missing argument value");
        std::string val=argv[++i];
        if(key=="--n")n=std::stoi(val);else if(key=="--steps")steps=std::stoi(val);
        else if(key=="--stride")stride=std::stoi(val);else if(key=="--threads")threads=std::stoi(val);
        else if(key=="--seed")seed=std::stoul(val);else if(key=="--density")density=std::stod(val);
        else if(key=="--temperature")temp=std::stod(val);else if(key=="--dt")dt=std::stod(val);
        else if(key=="--method")method=val;else if(key=="--out")out=val;else throw std::runtime_error("Unknown option: "+key);
    }
    if(n<2||steps<0||stride<1||threads<1||!std::isfinite(density)||density<=0||!std::isfinite(temp)||temp<=0||!std::isfinite(dt)||dt<=0)throw std::runtime_error("Invalid parameters");
    if(method!="grid"&&method!="sequential"&&method!="openmp")throw std::runtime_error("Unknown method");
    omp_set_num_threads(threads);Gas gas(n,density,temp,seed);
    std::filesystem::create_directories(out);
    std::ofstream traj(out+"/trajectory.csv"),obs(out+"/observables.csv"),meta(out+"/metadata.json");
    if(!traj||!obs||!meta)throw std::runtime_error("Cannot write output");
    traj<<std::setprecision(15)<<"step,id,x,y,vx,vy,xu,yu\n";obs<<std::setprecision(15)<<"step,time,kinetic,potential,total,temperature,px,py,pressure\n";
    gas.forces(method);double compute=0;
    for(int s=0;s<=steps;s++) {
        if(s%stride==0||s==steps) {
            double k=gas.kinetic();V p;for(auto w:gas.v){p.x+=w.x;p.y+=w.y;}
            if(!std::isfinite(k+gas.potential))throw std::runtime_error("Unstable integration: reduce dt or density");
            obs<<s<<','<<s*dt<<','<<k<<','<<gas.potential<<','<<k+gas.potential<<','<<k/(n-1)<<','<<p.x<<','<<p.y<<','<<(2*k+gas.virial)/(2*gas.L*gas.L)<<'\n';
            for(int i=0;i<n;i++)traj<<s<<','<<i<<','<<gas.r[i].x<<','<<gas.r[i].y<<','<<gas.v[i].x<<','<<gas.v[i].y<<','<<gas.unwrapped[i].x<<','<<gas.unwrapped[i].y<<'\n';
        }
        if(s<steps){auto start=std::chrono::steady_clock::now();gas.step(dt,method);compute+=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();}
    }
    meta<<std::setprecision(15)<<"{\"n\":"<<n<<",\"box\":"<<gas.L<<",\"density\":"<<density<<",\"initial_temperature\":"<<temp<<",\"dt\":"<<dt<<",\"steps\":"<<steps<<",\"seed\":"<<seed<<",\"method\":\""<<method<<"\",\"threads\":"<<threads<<",\"compute_seconds\":"<<compute<<"}\n";
    std::cout<<method<<": "<<compute<<" seconds; results: "<<out<<'\n';
}catch(const std::exception &e){std::cerr<<"Error: "<<e.what()<<'\n';return 1;}}
