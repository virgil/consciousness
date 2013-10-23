// Defines the unsigned double class
#include <assert.h>

class udouble
{
public:
    double dobb;
    uDouble(const double d) { dobb=d; }
    
    operator double() const {
        assert(dobb >= 0.0);
        return dobb;
    }
};

/*
void setTempo(uDouble dobb){ std::cout << dobb; }

int main()
{
    setTempo(3.14159);
    system("PAUSE");
    return 0;
} 
*/