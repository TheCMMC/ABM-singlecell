using namespace std;
#include <string>

extern "C" {
  void cfd_setup(char *);
  void cfd_query(double, double, double, double *, double *, double *);
};

int main(){
  char c[] = "foo";
  cfd_setup(c);
  double x=1,y=2,z=3,u,v,w;
  cfd_query(x,y,z,&u,&v,&w);
}
