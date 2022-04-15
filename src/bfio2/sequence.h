#include <cstdlib>
class Seq
{
    private:
        size_t start_index_, stop_index_, step_;
    public:
        inline Seq(size_t start, size_t  stop, size_t  step=1):start_index_(start), stop_index_(stop), step_(step){} 
        inline size_t Start()  const  {return start_index_;}
        inline size_t Stop()  const {return stop_index_;}
        inline size_t Step()  const {return step_;}
};