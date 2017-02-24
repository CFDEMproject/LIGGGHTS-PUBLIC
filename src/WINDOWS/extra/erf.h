# ifndef ERF_H
# define ERF_H

# if((_WIN32 || _WIN64) && _MSC_VER <= 1700)  //1800=Visual Studio 2013 has already erf and erfc


# ifdef __cplusplus
extern "C" {
# endif

double erf(double x);
double erfc(double x);

# ifdef __cplusplus
}
# endif

# endif

# endif
