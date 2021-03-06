#ifdef HAVE_CONFIG_H
 #include <config.hpp>
#endif

#include <chrono>
#include <fstream>
#include <iostream>
#include <sstream>
#include <map>
#include <cmath>
#include <vector>

using namespace std;

class sanfo_ostream
{
  /// buffered out
  ostream& out;
  
  /// buffer
  ostringstream buf;
  
public:
  
  void dump()
  {
    out<<buf.str();
    buf.clear();
    buf.str("");
  }
  
  /// constructor
  sanfo_ostream(ostream& out) : out(out)
  {
  }
  
  /// destructor
  ~sanfo_ostream()
  {
    dump();
  }
  
  /// print
  template <typename T>
  friend sanfo_ostream& operator<<(sanfo_ostream& os,const T& t)
  {
    os.buf<<t;
    
    return os;
  }
  
  /// print
  friend sanfo_ostream& operator<<(sanfo_ostream& os,const char t[])
  {
    int i=0;
    
    while(t[i]!='\0')
      {
	char c=t[i];
	os.buf<<c;
	if((c=='.' or c=='+') and os.buf.str().length()>100)
	  {
	    os.dump();
	    os.buf<<"\n";
	  }
	
	i++;
      }
    
    return os;
  }
};

sanfo_ostream fout(cout);

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

typedef int64_t Weight_t;

/// Separate blocks
constexpr char SEPARATOR[]=
  "/////////////////////////////////////////////////////////////////";

/// Identify an addend
using Addend=
  vector<int>;

/// Maps each addend to its weight
using WeightedAddend=
  map<Addend,Weight_t>;

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
sanfo_ostream& operator<<(sanfo_ostream& os,const Addend& a)
{
  ///Printing styles
  enum{FAPRILE,SUNPHO};
  
  /// Style
  constexpr bool printStyle=
    FAPRILE;
    
  for(int i=0;i<(int)a.size();i++)
    {
      if(i) os<<".";
      
      /// Val at pos i
      const int v=
	a[i];
      
      os<<(isM(v)?"m":"d")<<"["<<abs(v)<<"]";
    }
  
  // Append .d[0] if last entry is not d
  if constexpr(printStyle==FAPRILE)
    if(not isD(a.back()))
      os<<".d[0]";
  
  return os;
}

/// Returns an indented stream
sanfo_ostream& indent(sanfo_ostream& os,const int lev)
{
  for(int i=0;i<lev;i++)
    os<<" ";
  
  return os;
}

/// Print the result
sanfo_ostream& operator<<(sanfo_ostream& os,const WeightedAddend& x)
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
      const Weight_t aw=
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

/// Performs all manipulations and call iteratively
void simplifyAndAdd(WeightedAddend& res,Addend a,int pos=-1,const Weight_t w=1,int indentLev=0)
{
  stackSize++;
  stackSizeMax=max(stackSize,stackSizeMax);
  
  // If no pos was passed, starts from the last before end
  if(pos==-1)
    pos=a.size()-2;
  
  // Print input
  if constexpr(DEBUG)
    {
      indent(fout,indentLev)<<w<<"*"<<a<<"\n";
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
		  indent(fout,indentLev)<<"swapping"<<"\f";
		  indent(fout,indentLev)<<w<<"*"<<a<<"\f";
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
      
      // Write the current step
      if constexpr(DEBUG)
        indent(fout,indentLev)<<a<<"\n";
    }
  
  // Print some info on what is getting added
  if constexpr(DEBUG)
    {
      indent(fout,indentLev)<<"Adding "<<w<<" to "<<a<<" "<<res[a]<<"\n";
    cout<<"---"<<endl;
    }
  
  // Add the result
  res[a]+=w;
  
  // Decrease the stack size
  stackSize--;
}

/// Gets the n+1-th derivative from the n-th ine
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

/// Compute with a given n
void bench(const int n,const bool computeDirect=false)
{
  cout<<SEPARATOR<<" "<<n<<" "<<SEPARATOR<<endl;
  
  /////////////////////////////////////////////////////////////////
  
  /// Timer used for benchmarks
  double elapsedDirect=
    0;
  
  /// Result of the direct calculation
  WeightedAddend directRes;
  
  if(computeDirect)
    {
      directRes=benchmark(computeDirectly,elapsedDirect,n);
      
      cout<<directRes.size()<<" terms"<<endl;
      cout<<"Elapsed time to compute "<<n<<" directly: "<<elapsedDirect<<" s"<<endl;
      
      /////////////////////////////////////////////////////////////////
      cout<<SEPARATOR<<" "<<n<<" "<<SEPARATOR<<endl;
      /////////////////////////////////////////////////////////////////
    }
  
  /// Timer used for benchmarks
  double elapsedIter;
  
  /// Result of the iterative calculation
  WeightedAddend iterRes=
    benchmark(computeIterativetly,elapsedIter,n);
  
  cout<<iterRes.size()<<" terms"<<endl;
  cout<<"Elapsed time to compute "<<n<<" iteratively: "<<elapsedIter<<" s"<<endl;
  
  /////////////////////////////////////////////////////////////////
  
  /// Output file
  if constexpr(0)
    {
      ofstream rawIterResFile("ResIter"+to_string(n)+".txt");
      sanfo_ostream iterResFile(rawIterResFile);
      iterResFile<<iterRes<<"\n";
    }
  
  // Check and improvement report
  if(computeDirect)
    {
      fout<<iterRes-directRes<<"\n";
      
      cout<<"Improvement: "<<elapsedDirect/elapsedIter<<endl;
    }
}

int main(int narg,char **arg)
{
  // /// Degree of the derivative
  // const int n
  //   =30;
  
  // for(int i=6;i<=n;i++)
  //   bench(i);
  
  if(narg<2)
    CRASH("Use %s n",arg[0]);
  
  // Gets n as an argument
  int n;
  if(sscanf(arg[1],"%d",&n)!=1)
    CRASH("Error converting %s to n",arg[1]);
  
  fout<<computeIterativetly(n)<<"\n";
  
  return 0;
}
