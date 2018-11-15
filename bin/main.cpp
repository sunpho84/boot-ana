#ifdef HAVE_CONFIG_H
 #include <config.hpp>
#endif

#include <chrono>
#include <iostream>
#include <map>
#include <vector>

using namespace std;

/// Switch to enable debug
constexpr bool DEBUG=
  false;

/// Variadic print to a stream
///
/// Empty case
template <class S>
S& print(S& out)
{
  return out;
}

/// Variadic print to a stream
template <class S,class Head,class...Tail>
decltype(auto) print(S& out,const Head& head,Tail&&...tail)
{
  out<<head;
  
  if constexpr(sizeof...(tail)>0)
    {
      out<<" ";
      return print(out,forward<Tail>(tail)...);
    }
  else
    return out;
}

/// Macro to make the crash being more explicit
#define CRASH(...)                                                      \
  internalCrash(__LINE__,__FILE__,__PRETTY_FUNCTION__,__VA_ARGS__)

/// Crash with a detailed message
template <class...Args>
void internalCrash(const int line,const char *path,const char *funcName,const Args&...args)
{
  print(std::cerr,"ERROR in function ",funcName," at line ",line," of file ",path,": \"",args...,"\"\n");
  exit(1);
}

/// Get time
auto takeTime()
{
  return chrono::steady_clock::now();
}

/// Convert a time difference in second
double convTimeDiff(const decltype(takeTime()-takeTime())& diff)
{
  return chrono::duration<double>(diff).count();
}

/////////////////////////////////////////////////////////////////

/// Separate blocks
constexpr char SEPARATOR[]=
  "/////////////////////////////////////////////////////////////////";

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

/// True if is M or one of its derivative
bool isM(const int x)
{
  return x>=0;
}

/// True if is a derivative
bool isD(const int x)
{
  return x<0;
}

/// Print an addend
ostream& operator<<(ostream& os,const Addend& a)
{
  for(int i=0;i<(int)a.size();i++)
    {
      if(i) os<<".";
      
      /// Val at pos i
      const int v=
	a[i];
      
      os<<(isM(v)?"M":"D")<<"["<<abs(v)<<"]";
    }
  
  return os;
}

/// Returns an indented stream
ostream& indent(ostream& os,const int lev)
{
  for(int i=0;i<lev;i++)
    os<<" ";
  
  return os;
}

/// Print the result
ostream& operator<<(ostream& os,const WeightedAddend& x)
{
  bool first=
    true;
  
  // Print the result
  for(auto &[a,w] : x)
    {
      // Print the sign only if positive and not the first
      if((not first) and w>=0)
	 os<<"+";
      first=false;
      
      // Absolute value of the weight
      const int aw=
	abs(w);
      
      // Prints the coefficient only if not 1
      if(aw!=1)
	os<<abs(w)<<"*";
      
      os<<a;
    }
  
  return os;
}

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
	(i>>j)&0x1;
      
      ret[j]=
	b?M0:D1;
    }
  
  return ret;
}

/// Current stack size
int stackSize=
		    0;

/// Maximal dimension reached for the stack
int stackSizeMax=
		    0;

/// Performs all manipulaitions and call iteratively
void simplifyAndAdd(WeightedAddend& res,Addend a,int pos=-1,const int w=1,int indentLev=0)
{
  stackSize++;
  stackSizeMax=max(stackSize,stackSizeMax);
  
  // If no pos was passed, starts from the last before end
  if(pos==-1)
    pos=a.size()-2;
  
  // Print input
  if constexpr(DEBUG)
    {
      indent(cout,indentLev)<<w<<"*"<<a<<endl;
    }
  
  while(pos>=0 and a.size()!=1)
    {
      /// Reference to the value of the position
      int& v=
	a[pos];
      
      // if there is a D, move it right
      if(v==D1)
	{
	  /// Reference to the value of the next position
	  int& v_next=
	    a[pos+1];
	  
	  /// Check if next position is the last one
	  const bool nextIsLast=
	    pos+1==(int)a.size()-1;
	  
	  // If the next position is an M, swap and shift pos right if possible
	  if(isM(v_next))
	    {
	      /// Prepare the commutator
	      Addend c(a.size()-1);
	      
	      // Add until pos
	      for(int i=0;i<pos;i++) c[i]=a[i];
	      // Insert derivative of M incremented
	      c[pos]=a[pos+1]+1;
	      // Copy until the end
	      for(int i=pos+1;i<(int)c.size();i++) c[i]=a[i+1];
	      
	      // Call nested
	      simplifyAndAdd(res,c,pos,w,indentLev+1);
	      
	      swap(v_next,v);
	      if constexpr(DEBUG)
	        {
		  indent(cout,indentLev)<<"swapping"<<endl;
		  indent(cout,indentLev)<<w<<"*"<<a<<endl;
		}
	      
	      // If next position is the last one, do not move
	      if(not nextIsLast) pos++;
	    }
	  // If next is not an M, it's a D. Assuming is last one, sum it to current pos and remove last
	  else
	    {
	      if(not nextIsLast)
		CRASH("Next must be last");
	      
	      v+=v_next;
	      a.pop_back();
	    }
	}
      // If it's an M, ignore it
      else
	pos--;
      
      if constexpr(DEBUG)
        indent(cout,indentLev)<<a<<endl;
    }
  
  // Print some info on what is getting added
  if constexpr(DEBUG)
    {
      indent(cout,indentLev)<<"Adding "<<w<<" to "<<a<<" "<<res[a]<<endl;
    cout<<"---"<<endl;
    }
  
  // Add the result
  res[a]+=w;
  
  // Decrease the stack size
  stackSize--;
}

/// Gets the n+1 from n
WeightedAddend differentiate(const WeightedAddend& in)
{
  /// Result
  WeightedAddend out;
  
  for(const auto& [a,w] : in)
    {
      for(const int ins : {M0,D1})
	{
	  /// Term with n
	  Addend n(1);
	  n.reserve(a.size()+1);
	  n[0]=ins;
	  n.insert(n.begin()+1,a.begin(),a.end());
	  
	  simplifyAndAdd(out,n,0,w);
	}
    }
  
  return out;
}

/// Differentiate n times
WeightedAddend differentiate(const WeightedAddend& in,const int n)
{
  /// Result
  WeightedAddend out=
    in;
  
  for(int i=0;i<n;i++)
   out=differentiate(out);
  
 return out;
}

/// Compute the derivative iteratively
WeightedAddend computeIterativetly(const int n)
{
  /// Null polynomial
  WeightedAddend zero{{Addend{},1}};
  
  return differentiate(zero,n);
}

/// Compute the derivative diectly
WeightedAddend computeDirectly(const int n)
{
  /// Result
  WeightedAddend res;
  
  /// Number of term
  const int nAddend=
    1<<n;
  
  // Loop on each addend
  for(int i=0;i<nAddend;i++)
    simplifyAndAdd(res,getTerm(i,n));
  
  return res;
}

/// Compute printing the time
template <typename F,
	  typename...Args>
decltype(auto) benchmark(F f,double &elaps,Args&&...args)
{
  auto beg=
    takeTime();
  
  auto res=
    f(forward<Args>(args)...);
  
  auto end=
    takeTime();
  
  elaps=convTimeDiff(end-beg);
  
  return res;
}

/// Removes all zero elements of in
WeightedAddend prune(const WeightedAddend& in)
{
  WeightedAddend res;
  
  for(auto& [a,w] : in)
    if(w!=0)
      res[a]=w;
  
  return res;
}

/// Takes the difference between a and b
WeightedAddend operator-(const WeightedAddend& in1,const WeightedAddend& in2)
{
  WeightedAddend res=
    in1;
  
  for(auto& [a,w] : in2)
    res[a]-=w;
  
  return prune(res);
}

int main(int narg,char **arg)
{
  /// Degree of the derivative
  const int n
    =3;
  
  cout<<SEPARATOR<<" "<<n<<" "<<SEPARATOR<<endl;
  
  /// Timer used for benchmarks
  double elapsed;
  
  /////////////////////////////////////////////////////////////////
  
  WeightedAddend directRes=
    benchmark(computeDirectly,elapsed,n);
  
  cout<<directRes<<endl;
  cout<<"Elapsed time to compute directly: "<<elapsed<<" s"<<endl;
  
  /////////////////////////////////////////////////////////////////
  cout<<SEPARATOR<<" "<<n<<" "<<SEPARATOR<<endl;
  /////////////////////////////////////////////////////////////////
  
  WeightedAddend iterRes=
    benchmark(computeIterativetly,elapsed,n);
  
  cout<<iterRes<<endl;
  cout<<"Elapsed time to compute iteratively: "<<elapsed<<" s"<<endl;
  
  cout<<iterRes-directRes<<endl;
  
  return 0;
}
