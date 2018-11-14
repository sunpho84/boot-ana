#ifdef HAVE_CONFIG_H
 #include <config.hpp>
#endif

#include <map>
#include <vector>

using namespace std;

/// Identify an addend
using Addend=
  vector<int>;

/// Maps each addend to its weight
using WeightedAddend=
  map<Addend,int>;

/// Representation of M
constexpr int M0=
  0;

/// Representation of D
constexpr int D1=
  -1;

/// Associate the i-th of the n^2 terms in the (D+m)^n expansion
///
/// Decompose the bit representation, and maps 0 to D, 1 to m
Addend getTerm(const int i,const int n)
{
  /// Returned value
  Addend ret(n);
  
  for(int j=0;j<n;j++)
    {
      /// j-th bit of i
      const bool b=
	(i>>j)%0x1;
      
      ret[j]=
	b?M0:D1;
    }
  
  return ret;
}

void simplfyAndAdd()
{
}

int main(int narg,char **arg)
{
  const int n=
    10;
  
  
  return 0;
}
