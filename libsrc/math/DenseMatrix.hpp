#pragma once
#include <vector>
namespace prl
{
  using real_t = double;
  using Vector = std::vector<real_t>;
  
  struct DenseMatrix
  {
    //ctor
    DenseMatrix()
      {
        //nothing to do..
      }
    
    //nxn matrix.
    DenseMatrix(int n)
    {
      Init(n);
    }
    // nxm matrix.
    DenseMatrix(int n, int m)
    {
      Init(n,m);
    }

    void Init(int n)
    {
      ncols_ = n;
      nrows_ = n;
      nents_ = n*n;

      A_.clear();
      A_.resize(nents_);
      std::fill(A_.begin(),A_.end(),0.0);
    }

    // nxm matrix.
    void Init(int n, int m)
    {
      nrows_=n;
      ncols_=m;
      is_square_ = false;
      nents_ = n*m;

      A_.clear();
      A_.resize(nents_);
      std::fill(A_.begin(),A_.end(),0.0);
    }

    void SetSize(int n)
    {
      Init(n);
    }
    
    int Width() const
    {
      return ncols_;
    }
    int Height() const
    {
      return nrows_;
    }
    
    const real_t & operator()(const int n, const int m) const
    {
      return A_[n*ncols_+m];
    }
    real_t & operator()(const int n, const int m) 
    {
      return A_[n*ncols_+m];
    }
    
    DenseMatrix & operator = (const real_t val) 
    {
      std::fill(A_.begin(),A_.end(),val);
      return *this;
      
    }
    
    bool is_square_ = true;
    int ncols_;
    int nrows_;
    int nents_;
    std::vector<real_t> A_;
    
      
  };//struct DenseMatrix
// helper functions
  // y=Ax;
  inline void Matvec(const DenseMatrix &A, const Vector &x, Vector&y) 
  {
    int width=A.Width();
    int height=A.Height();
    for(int i=0;i<height;++i)
    {
      y[i]=0;
      for(int j=0;j<width;++j)
      {
        y[i] += A(i,j)*x[j];
      }//j
    }//i
  }//Matvec

  inline void Scale(const real_t factor, Vector&x) 
    {
      int size=x.size();
      for(int i=0;i<size;++i)
        x[i] *=factor;
    }
  inline void Scale(const real_t factor, DenseMatrix &A) 
    {
      int width=A.Width();
      int height= A.Height();
      
      for(int i=0;i<height;++i)
        for(int j=0;j<width;++j)
          A(i,j) *=factor;
    }//Scale

  inline void Add(const real_t factor, const Vector &vec, Vector& result) 
    {
      int size=result.size();
      
      for( int i=0;i<size;++i)
        result[i] += factor*vec[i];
      
    }
  

}//namespace prl
