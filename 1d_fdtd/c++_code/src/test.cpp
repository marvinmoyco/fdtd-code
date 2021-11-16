#include <xtensor-fftw/basic.hpp>   // rfft, irfft
#include <xtensor-fftw/helper.hpp>  // rfftscale
#include <xtensor/xarray.hpp>
#include <xtensor/xbuilder.hpp>     // xt::arange
#include <xtensor/xmath.hpp>        // xt::sin, cos
#include <complex>
#include <xtensor/xio.hpp>

int main()
{
    // generate a sinusoid field
    double dx = M_PI / 100;
    xt::xarray<double> x = xt::arange(0., 2 * M_PI, dx);
    xt::xarray<double> sin_x = xt::sin(x);

    // transform to Fourier space
    auto sin_fs = xt::fftw::rfft(sin_x);
    std::cout << "Shape of x: (" << x.shape()[0] << "," << x.shape()[1] << ")" << std::endl;
    std::cout << "Shape of sin: (" << sin_x.shape()[0] << "," << sin_x.shape()[1] << ")" << std::endl;
    std::cout << "Shape of sin_fs: (" << sin_fs.shape()[0] << "," << sin_fs.shape()[1] << ")" << std::endl;
    // multiply by i*k
    std::complex<double> i {0, 1};
    auto k = xt::fftw::rfftscale<double>(sin_x.shape()[0], dx);
    xt::xarray<std::complex<double>> sin_derivative_fs = xt::eval(i * k * sin_fs);

    // transform back to normal space
    auto sin_derivative = xt::fftw::irfft(sin_derivative_fs);
    std::cout << "Shape of sin_derivative: (" << sin_derivative.shape()[0] << "," << sin_derivative.shape()[1] << ")" << std::endl;

    return 0;
}