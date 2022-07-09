#ifndef FORMATTER_HPP
#define FORMATTER_HPP

#include <sstream>
#include <vector>

namespace prl {

class formatter_t
{
public:
    formatter_t() {}
    ~formatter_t() {}

    template <typename Type>
    formatter_t & operator << (const Type & value)
    {
        stream_ << value;
        return *this;
    }
    
    template<typename Type>
    formatter_t & operator << (const std::vector<Type> & values)
    {
      if (!values.empty()) {
        for (unsigned i=0; i<values.size()-1; ++i) stream_ << values[i] << ", ";
        stream_ << values.back();
      }
      return *this;
    }

    std::string str() const         { return stream_.str(); }
    operator std::string () const   { return stream_.str(); }

    enum ConvertToString 
    {
        to_str
    };
    std::string operator >> (ConvertToString) { return stream_.str(); }

private:
    std::stringstream stream_;

    formatter_t(const formatter_t &);
    formatter_t & operator = (formatter_t &);
};

} // namespace

#endif // CONTRA_FORMATTER_HPP
