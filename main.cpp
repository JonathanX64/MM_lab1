#include <iostream>
#include <thread>
#include <future>
#include <unistd.h>
#include "bitmap_image.hpp"
#include <fcntl.h>
#include <math.h>

bitmap_image leave_single_color(const bitmap_image &source, char color);

bitmap_image mirror(const bitmap_image &source, char direction);

bitmap_image generate_rgb_back(const bitmap_image &y, const bitmap_image &cb, const bitmap_image &cr);

bitmap_image rotate_ccw(const bitmap_image &source);

bitmap_image resizeBilinear(const bitmap_image &source, unsigned int w2, unsigned int h2);

bitmap_image decimation(const bitmap_image &source);

void histograms(const bitmap_image &image);

void rotation_cycle(const bitmap_image &image);

bitmap_image histo(const bitmap_image &source, char channel);

bitmap_image histo(std::vector<std::vector<int>> &source);

void ycbcr_histos(const bitmap_image &y, const bitmap_image &cb, const bitmap_image &cr);

double count_expected_value(const bitmap_image &source, char component);

double count_correlation(const bitmap_image &source, char component1, char component2);

double count_standard_deviation(const bitmap_image &source, char component);

double enthropy(const bitmap_image &rgb, char channel);

bitmap_image merge_ycbcr(const bitmap_image &y, const bitmap_image &cb, const bitmap_image &cr);

std::vector<std::vector<std::vector<int>>> DPCM(const bitmap_image &source, char channel);

std::vector<double> dpcm_enthropies(const std::vector<std::vector<std::vector<int>>> &data, const bitmap_image &source);

std::vector<double> DPCM_driver(const bitmap_image &source, char channel);

std::vector<double> DPCM_driver(const bitmap_image &source, char channel) {
    std::vector<std::vector<std::vector<int>>> data = DPCM(source, channel);
    std::vector<double> entropies = dpcm_enthropies(data, source);
    return entropies;
}

void dpcm_histos(const bitmap_image &source, const bitmap_image &merged_ycbcr) {
    std::vector<std::vector<std::vector<int>>> data1 = DPCM(source, 'r');
    std::vector<std::vector<std::vector<int>>> data2 = DPCM(source, 'g');
    std::vector<std::vector<std::vector<int>>> data3 = DPCM(source, 'b');
    std::vector<std::vector<std::vector<int>>> data4 = DPCM(merged_ycbcr, 'r');
    std::vector<std::vector<std::vector<int>>> data5 = DPCM(merged_ycbcr, 'g');
    std::vector<std::vector<std::vector<int>>> data6 = DPCM(merged_ycbcr, 'b');

    histo(data1[0]).save_image("images/15. R1.bmp");
    histo(data1[1]).save_image("images/15. R2.bmp");
    histo(data1[2]).save_image("images/15. R3.bmp");
    histo(data1[3]).save_image("images/15. R4.bmp");

    histo(data2[0]).save_image("images/15. G1.bmp");
    histo(data2[1]).save_image("images/15. G2.bmp");
    histo(data2[2]).save_image("images/15. G3.bmp");
    histo(data2[3]).save_image("images/15. G4.bmp");

    histo(data3[0]).save_image("images/15. B1.bmp");
    histo(data3[1]).save_image("images/15. B2.bmp");
    histo(data3[2]).save_image("images/15. B3.bmp");
    histo(data3[3]).save_image("images/15. B4.bmp");

    histo(data4[0]).save_image("images/15. Y1.bmp");
    histo(data4[1]).save_image("images/15. Y2.bmp");
    histo(data4[2]).save_image("images/15. Y3.bmp");
    histo(data4[3]).save_image("images/15. Y4.bmp");

    histo(data5[0]).save_image("images/15. Cb1.bmp");
    histo(data5[1]).save_image("images/15. Cb2.bmp");
    histo(data5[2]).save_image("images/15. Cb3.bmp");
    histo(data5[3]).save_image("images/15. Cb4.bmp");

    histo(data6[0]).save_image("images/15. Cr1.bmp");
    histo(data6[1]).save_image("images/15. Cr2.bmp");
    histo(data6[2]).save_image("images/15. Cr3.bmp");
    histo(data6[3]).save_image("images/15. Cr4.bmp");
}

void DPCM_cycle(const bitmap_image &source, const bitmap_image &merged_ycbcr) {
    std::future<std::vector<double>> i1 = std::async(&DPCM_driver, std::ref(source), 'r');
    std::future<std::vector<double>> i2 = std::async(&DPCM_driver, std::ref(source), 'g');
    std::future<std::vector<double>> i3 = std::async(&DPCM_driver, std::ref(source), 'b');
    std::future<std::vector<double>> i4 = std::async(&DPCM_driver, std::ref(merged_ycbcr), 'r');
    std::future<std::vector<double>> i5 = std::async(&DPCM_driver, std::ref(merged_ycbcr), 'g');
    std::future<std::vector<double>> i6 = std::async(&DPCM_driver, std::ref(merged_ycbcr), 'b');

    std::vector<double> result1 = i1.get();
    std::vector<double> result2 = i2.get();
    std::vector<double> result3 = i3.get();
    std::vector<double> result4 = i4.get();
    std::vector<double> result5 = i5.get();
    std::vector<double> result6 = i6.get();

    std::cout << "\n16. Enthropy of DPCM array of R channel" << std::endl;
    for (int i = 0; i < 4; ++i) {
        std::cout << i + 1 << " way: ";
        std::cout << "H(D) = " << result1[i] << std::endl;
    }

    std::cout << "Enthropy of DPCM array of G channel" << std::endl;
    for (int i = 0; i < 4; ++i) {
        std::cout << i + 1 << " way: ";
        std::cout << "H(D) = " << result2[i] << std::endl;
    }

    std::cout << "Enthropy of DPCM array of B channel" << std::endl;
    for (int i = 0; i < 4; ++i) {
        std::cout << i + 1 << " way: ";
        std::cout << "H(D) = " << result3[i] << std::endl;
    }

    std::cout << "Enthropy of DPCM array of Y channel" << std::endl;
    for (int i = 0; i < 4; ++i) {
        std::cout << i + 1 << " way: ";
        std::cout << "H(D) = " << result4[i] << std::endl;
    }

    std::cout << "Enthropy of DPCM array of Cb channel" << std::endl;
    for (int i = 0; i < 4; ++i) {
        std::cout << i + 1 << " way: ";
        std::cout << "H(D) = " << result5[i] << std::endl;
    }

    std::cout << "Enthropy of DPCM array of Cr channel" << std::endl;
    for (int i = 0; i < 4; ++i) {
        std::cout << i + 1 << " way: ";
        std::cout << "H(D) = " << result6[i] << std::endl;
    }
}

void call_from_thread(const bitmap_image &image, bitmap_image &result, char color, const std::string &path) {
    result = leave_single_color(image, color);
    result.save_image(path);
}

void ycbcr_histos(const bitmap_image &y, const bitmap_image &cb, const bitmap_image &cr) {
    histo(y, 'r').save_image("images/13. Y hist.bmp");
    histo(cb, 'r').save_image("images/13. Cb hist.bmp");
    histo(cr, 'r').save_image("images/13. Cr hist.bmp");
}

void correlations(const bitmap_image &image, const bitmap_image &y, const bitmap_image &cb, const bitmap_image &cr) {
    std::future<double> i1 = std::async(&count_correlation, image, 'r', 'g');
    std::future<double> i2 = std::async(&count_correlation, image, 'r', 'b');
    std::future<double> i3 = std::async(&count_correlation, image, 'b', 'g');

    std::cout << "4. Correlation coefficient between R n G: " << i1.get() << std::endl;
    std::cout << "Correlation coefficient between R n B: " << i2.get() << std::endl;
    std::cout << "Correlation coefficient between B n G: " << i3.get() << std::endl;

    bitmap_image tmp = merge_ycbcr(y, cb, cr);
    std::future<double> i4 = std::async(&count_correlation, tmp, 'r', 'g');
    std::future<double> i5 = std::async(&count_correlation, tmp, 'r', 'b');
    std::future<double> i6 = std::async(&count_correlation, tmp, 'b', 'g');
    std::cout << "\n5. Correlation coefficient between Y n Cb: " << i4.get() << std::endl;
    std::cout << "Correlation coefficient between Cb n Cr: " << i5.get() << std::endl;
    std::cout << "Correlation coefficient between Cr n Y: " << i6.get() << std::endl;
}

void enthropies(const bitmap_image &image, const bitmap_image &y, const bitmap_image &cb, const bitmap_image &cr) {
    std::future<double> i1 = std::async(&enthropy, image, 'r');
    std::future<double> i2 = std::async(&enthropy, image, 'g');
    std::future<double> i3 = std::async(&enthropy, image, 'b');

    std::cout << "\n13. Enthropy of R channel: " << i1.get() << std::endl;
    std::cout << "Enthropy of G channel: " << i2.get() << std::endl;
    std::cout << "Enthropy of B channel: " << i3.get() << std::endl;

    std::future<double> i4 = std::async(&enthropy, y, 'r');
    std::future<double> i5 = std::async(&enthropy, cb, 'r');
    std::future<double> i6 = std::async(&enthropy, cr, 'r');
    std::cout << "Enthropy of Y channel: " << i4.get() << std::endl;
    std::cout << "Enthropy of Cb channel: " << i5.get() << std::endl;
    std::cout << "Enthropy of Cr channel: " << i6.get() << std::endl;
}

int main(int argc, char **argv) {

    std::string path;
    if (argc == 2) {
        path = argv[1];
    } else {
        path = "images/a.bmp";
    }

    bitmap_image image(path);
    std::thread t8(histograms, std::ref(image));

    bitmap_image y, cb, cr;
    std::thread t4(call_from_thread, std::ref(image), std::ref(y), 'y', "images/6. Y only.bmp");
    std::thread t5(call_from_thread, std::ref(image), std::ref(cb), 'B', "images/6. Cb only.bmp");
    std::thread t6(call_from_thread, std::ref(image), std::ref(cr), 'R', "images/6. Cr only.bmp");

    t4.join();
    t5.join();
    t6.join();
    dpcm_histos(image, merge_ycbcr(y, cb, cr)); //15
    std::thread t10(correlations, std::ref(image), std::ref(y), std::ref(cb), std::ref(cr));

    bitmap_image r, g, b;
    std::thread t1(call_from_thread, std::ref(image), std::ref(r), 'r', "images/3. Red only.bmp");
    std::thread t2(call_from_thread, std::ref(image), std::ref(g), 'g', "images/3. Green only.bmp");
    std::thread t3(call_from_thread, std::ref(image), std::ref(b), 'b', "images/3. Blue only.bmp");

    t1.join();
    t2.join();
    t3.join();

    std::thread t7(rotation_cycle, std::ref(image));
    t8.join();

    std::thread t9(ycbcr_histos, std::ref(y), std::ref(cb), std::ref(cr));

    bitmap_image ycbcr_restored = generate_rgb_back(y, cb, cr);
    ycbcr_restored.save_image("images/7. Restored after RGB-YCbCr-RGB.bmp");

    t10.join();
    std::cout << "\n7. PSNR of r channel after RGB->YCBCR->RGB conversion: " << ycbcr_restored.channel_psnr(image, 'r')
              << std::endl;
    std::cout << "PSNR of g channel after RGB->YCBCR->RGB conversion: " << ycbcr_restored.channel_psnr(image, 'g')
              << std::endl;
    std::cout << "PSNR of b channel after RGB->YCBCR->RGB conversion: " << ycbcr_restored.channel_psnr(image, 'b')
              << std::endl;

    bitmap_image x_mirrored = mirror(image, 'w'), y_mirrored = mirror(image, 'h');
    x_mirrored.save_image("images/1. Mirrored by x axis.bmp");
    y_mirrored.save_image("images/1. Mirrored by y axis.bmp");

    t7.join();

    bitmap_image deci_cr = decimation(cr);
    deci_cr.save_image("images/8a. decimated_cr.bmp");

    bitmap_image deci_cr_2 = resizeBilinear(cr, cr.width() / 2, cr.height() / 2);
    deci_cr_2.save_image("images/8b. decimated_cr_2.bmp");

    bitmap_image deci_cb = decimation(cb);
    deci_cb.save_image("images/8a. decimated_cb.bmp");

    bitmap_image deci_cb_2 = resizeBilinear(cb, cb.width() / 2, cb.height() / 2);
    deci_cb_2.save_image("images/8b. decimated_cb_2.bmp");

    t9.join();

    bitmap_image cb_restored = resizeBilinear(deci_cb_2, cb.width(), cb.height());
    bitmap_image cr_restored = resizeBilinear(deci_cr_2, cr.width(), cr.height());
    bitmap_image after_restoraion = generate_rgb_back(y, cb_restored, cr_restored);
    after_restoraion.save_image("images/9. Restored after RGB-YCbCr-50% scale-RGB.bmp");

    std::cout << "\n10. PSNR of r channel after RGB->YCBCR->50% scale->RGB conversion: "
              << after_restoraion.channel_psnr(image, 'r')
              << std::endl;
    std::cout << "PSNR of g channel after RGB->YCBCR->50% scale->RGB conversion: "
              << after_restoraion.channel_psnr(image, 'g')
              << std::endl;
    std::cout << "PSNR of b channel after RGB->YCBCR->50% scale->RGB conversion: "
              << after_restoraion.channel_psnr(image, 'b')
              << std::endl;

    bitmap_image deci_cr_3 = resizeBilinear(cr, cr.width() / 4, cr.height() / 4);
    deci_cr_3.save_image("images/11. decimated_cr_3.bmp");

    bitmap_image deci_cb_3 = resizeBilinear(cb, cb.width() / 4, cb.height() / 4);
    deci_cb_3.save_image("images/11. decimated_cb_3.bmp");

    bitmap_image cb_restored_2 = resizeBilinear(deci_cb_3, cb.width(), cb.height());
    bitmap_image cr_restored_2 = resizeBilinear(deci_cr_3, cr.width(), cr.height());
    bitmap_image after_another_restoraion = generate_rgb_back(y, cb_restored_2, cr_restored_2);
    after_restoraion.save_image("images/11. Restored after RGB-YCbCr-25% scale-RGB.bmp");

    std::cout << "\n11. PSNR of r channel after RGB->YCBCR->25% scale->RGB conversion: "
              << after_another_restoraion.channel_psnr(image, 'r')
              << std::endl;
    std::cout << "PSNR of g channel after RGB->YCBCR->25% scale->RGB conversion: "
              << after_another_restoraion.channel_psnr(image, 'g')
              << std::endl;
    std::cout << "PSNR of b channel after RGB->YCBCR->25% scale->RGB conversion: "
              << after_another_restoraion.channel_psnr(image, 'b') << std::endl;
    enthropies(image, y, cb, cr);

    DPCM_cycle(image, merge_ycbcr(y, cb, cr)); //16
    return 0;
}

//TODO .4.2

double enthropy(const bitmap_image &rgb, char channel) {
    std::vector<double> result(256);
    for (int i = 0; i < 256; i++) {
        result[i] = 0;
    }

    for (uint16_t i = 0; i < rgb.width(); i++) {
        for (uint16_t j = 0; j < rgb.height(); j++) {
            switch (channel) {
                case 'r':
                    result[rgb.get_pixel(i, j).red]++;
                    break;
                case 'g':
                    result[rgb.get_pixel(i, j).green]++;
                    break;
                case 'b':
                    result[rgb.get_pixel(i, j).blue]++;
                    break;
                default:
                    std::cout << "default case reached! Something went wrong..." << std::endl;
                    break;
            }
        }
    }

    for (int j = 0; j < 256; j++) {
        result[j] /= rgb.height() * rgb.width();
    }

    double H = 0;
    for (int i = 0; i < 256; i++) {
        if (result[i] != 0) {
            H -= result[i] * log2(result[i]);
        }
    }
    return H;
}

double count_correlation(const bitmap_image &source, char component1, char component2) { //коэф корреляции
    double r = 0;
    double tmp = 0;
    unsigned char A;
    unsigned char B;
    double M1 = 0;
    double M2 = 0;

    switch (component1) {
        case 'r':
            M1 = count_expected_value(source, 'r');
            break;
        case 'b':
            M1 = count_expected_value(source, 'b');
            break;
        case 'g':
            M1 = count_expected_value(source, 'g');
            break;
        default:
            M1 = 0;
            std::cout << "default case reached! Something went wrong..." << std::endl;
            break;
    }

    switch (component2) {
        case 'r':
            M2 = count_expected_value(source, 'r');
            break;
        case 'b':
            M2 = count_expected_value(source, 'b');
            break;
        case 'g':
            M2 = count_expected_value(source, 'g');
            break;
        default:
            M2 = 0;
            std::cout << "default case reached! Something went wrong..." << std::endl;
            break;
    }

    double sigma = count_standard_deviation(source, component1) * count_standard_deviation(source, component2);
    for (int i = 0; i < source.height(); i++) {
        for (int j = 0; j < source.width(); j++) {
            switch (component1) {
                case 'r':
                    A = source.get_pixel(j, i).red;
                    break;
                case 'b':
                    A = source.get_pixel(j, i).blue;
                    break;
                case 'g':
                    A = source.get_pixel(j, i).green;
                    break;
                default:
                    A = 0;
                    std::cout << "default case reached! Something went wrong..." << std::endl;
                    break;
            }

            switch (component2) {
                case 'r':
                    B = source.get_pixel(j, i).red;
                    break;
                case 'b':
                    B = source.get_pixel(j, i).blue;
                    break;
                case 'g':
                    B = source.get_pixel(j, i).green;
                    break;
                default:
                    B = 0;
                    std::cout << "default case reached! Something went wrong..." << std::endl;
                    break;
            }

            tmp += (int(A) - M1) * (int(B) - M2);
        }
    }
    tmp /= source.width() * source.height();
    r = tmp / sigma;
    return r;
}

double count_standard_deviation(const bitmap_image &source, char component) { //Средн кв отклонение

    double sigma = 0;
    double M_red = count_expected_value(source, 'r');
    double M_blue = count_expected_value(source, 'b');
    double M_green = count_expected_value(source, 'g');

    for (int i = 0; i < source.width(); i++) {
        for (int j = 0; j < source.height(); j++) {
            switch (component) {
                case 'r':
                    sigma += pow(source.get_pixel(i, j).red - M_red, 2);
                    break;
                case 'g':
                    sigma += pow(source.get_pixel(i, j).green - M_green, 2);
                    break;
                case 'b':
                    sigma += pow(source.get_pixel(i, j).blue - M_blue, 2);
                    break;
                default:
                    sigma = 0;
                    std::cout << "default case reached! Something went wrong..." << std::endl;
                    break;
            }
        }
    }
    sigma = sqrt(sigma / (source.height() * source.width() - 1));
    return sigma;
}

double count_expected_value(const bitmap_image &source, char component) { //Мат ожидание
    double M = 0;
    for (int i = 0; i < source.width(); i++) {
        for (int j = 0; j < source.height(); j++) {
            switch (component) {
                case 'r':
                    M += source.get_pixel(i, j).red;
                    break;
                case 'g':
                    M += source.get_pixel(i, j).green;
                    break;
                case 'b':
                    M += source.get_pixel(i, j).blue;
                    break;
                default:
                    M += 0;
                    std::cout << "default case reached! Something went wrong..." << std::endl;
                    break;
            }

        }
    }
    M /= (source.width() * source.height());
    return M;
}

bitmap_image histo(const bitmap_image &source, char channel) {
    bitmap_image res(256, 256);
    double max = 0;
    uint32_t pixels[256];

    unsigned int current;
    for (uint32_t &pixel : pixels) {
        pixel = 0;
    }
    for (unsigned int i = 0; i < source.width(); ++i) {
        for (unsigned int j = 0; j < source.height(); ++j) {
            switch (channel) {
                case 'r':
                    current = source.get_pixel(i, j).red;
                    break;
                case 'g':
                    current = source.get_pixel(i, j).green;
                    break;
                case 'b':
                    current = source.get_pixel(i, j).blue;
                    break;
                case 'a':
                    current = static_cast<unsigned int>(
                            (source.get_pixel(i, j).red +
                             source.get_pixel(i, j).green +
                             source.get_pixel(i, j).blue) /
                            3.0);
                    break;
                default:
                    current = 0;
                    std::cout << "default case reached! Something went wrong..." << std::endl;
                    break;
            }

            pixels[current]++;
            if (pixels[current] > max) max = pixels[current];
        }
    }

    for (unsigned int i = 0; i < 256; ++i) {
        for (unsigned int j = 0; j < ((pixels[i] / max) * 256.0); ++j) {
            res.set_pixel(i, j, 126, 232, 255);
        }
    }

    return mirror(res, 'w');
}

bitmap_image histo(std::vector<std::vector<int>> &source) {
    bitmap_image res(256, 256);
    double max = 0;
    uint32_t pixels[256];

    unsigned int current;
    for (uint32_t &pixel : pixels) {
        pixel = 0;
    }

    for (unsigned int i = 0; i < source.size(); ++i) {
        for (unsigned int j = 0; j < source[0].size(); ++j) {
            current = source[i][j] % 255;
            if (current > 256) break;
            pixels[current]++;
            if (pixels[current] > max) max = pixels[current];
        }
    }
    //std::cout << "first cycle done" << std::endl;

    for (unsigned int i = 0; i < 256; ++i) {
        for (unsigned int j = 0; j < ((pixels[i] / max) * 255.0); ++j) {
            res.set_pixel(i, j, 176, 184, 65);
        }
    }

    //std::cout << "second cycle done" << std::endl;

    return mirror(res, 'w');
}

bitmap_image decimation(const bitmap_image &source) {
    bitmap_image res(source.width() / 2, source.height() / 2);
    for (unsigned int i = 0; i < res.width(); ++i) {
        for (unsigned int j = 0; j < res.height(); ++j) {
            res.set_pixel(i, j, source.get_pixel(i * 2, j * 2));
        }
    }

    return res;
}

bitmap_image merge_ycbcr(const bitmap_image &y, const bitmap_image &cb, const bitmap_image &cr) {
    bitmap_image res(y.width(), y.height());
    for (unsigned int i = 0; i < res.width(); ++i) {
        for (unsigned int j = 0; j < res.height(); ++j) {
            res.set_pixel(i, j,
                          y.get_pixel(i, j).red,
                          cb.get_pixel(i, j).red, //green
                          cr.get_pixel(i, j).red  //blue
            );
        }
    }
    return res;
}

void rotation_cycle(const bitmap_image &image) {
    bitmap_image result = rotate_ccw(image);
    result.save_image("images/1. Rotated_270.bmp");
    result = rotate_ccw(result);
    result.save_image("images/1. Rotated_180.bmp");
    result = rotate_ccw(result);
    result.save_image("images/1. Rotated_90.bmp");
}

void histograms(const bitmap_image &image) {
    histo(image, 'r').save_image("images/13. R hist.bmp");
    histo(image, 'g').save_image("images/13. G hist.bmp");
    histo(image, 'b').save_image("images/13. B hist.bmp");
    histo(image, 'a').save_image("images/13. RGB hist.bmp");
}

bitmap_image resizeBilinear(const bitmap_image &source, const unsigned int w2,
                            const unsigned int h2) {
    const unsigned int w = source.width();
    const unsigned int h = source.height();

    uint32_t length = source.width() * source.height() * 3;
    int *pixels = new int[source.width() * source.height()];
    uint32_t tmp = 0;
    for (uint32_t i = 0, k = 0; i < length; i += 3, k++) {
        tmp += 0xFF;
        tmp = tmp << 8u;

        tmp += source.data()[i];
        tmp = tmp << 8u;

        tmp += source.data()[i + 1];
        tmp = tmp << 8u;

        tmp += source.data()[i + 2];

        pixels[k] = tmp;
    }

    bitmap_image res(w2, h2);
    unsigned int a, b, c, d, x, y, index;
    float x_ratio = ((float) (w - 1)) / w2;
    float y_ratio = ((float) (h - 1)) / h2;
    float x_diff, y_diff;
    rgb_t pixel;
    for (uint32_t i = 0; i < h2; i++) {
        for (uint32_t j = 0; j < w2; j++) {
            x = (int) (x_ratio * j);
            y = (int) (y_ratio * i);
            x_diff = (x_ratio * j) - x;
            y_diff = (y_ratio * i) - y;
            index = (y * w + x);
            a = pixels[index];
            b = pixels[index + 1];
            c = pixels[index + w];
            d = pixels[index + w + 1];

            // Yb = Ab(1-w)(1-h) + Bb(w)(1-h) + Cb(h)(1-w) + Db(wh)
            pixel.blue = static_cast<unsigned char>((a & 0xffu) * (1 - x_diff) * (1 - y_diff) +
                                                    (b & 0xffu) * (x_diff) * (1 - y_diff) +
                                                    (c & 0xffu) * (y_diff) * (1 - x_diff) +
                                                    (d & 0xffu) * (x_diff * y_diff));

            // Yg = Ag(1-w)(1-h) + Bg(w)(1-h) + Cg(h)(1-w) + Dg(wh)
            pixel.green = static_cast<unsigned char>(((a >> 8u) & 0xff) * (1 - x_diff) * (1 - y_diff) +
                                                     ((b >> 8u) & 0xff) * (x_diff) * (1 - y_diff) +
                                                     ((c >> 8u) & 0xff) * (y_diff) * (1 - x_diff) +
                                                     ((d >> 8u) & 0xff) * (x_diff * y_diff));

            // Yr = Ar(1-w)(1-h) + Br(w)(1-h) + Cr(h)(1-w) + Dr(wh)
            pixel.red = static_cast<unsigned char>(((a >> 16u) & 0xff) * (1 - x_diff) * (1 - y_diff) +
                                                   ((b >> 16u) & 0xff) * (x_diff) * (1 - y_diff) +
                                                   ((c >> 16u) & 0xff) * (y_diff) * (1 - x_diff) +
                                                   ((d >> 16u) & 0xff) * (x_diff * y_diff));

            res.set_pixel(j, i, pixel);
        }
    }
    return res;
}

bitmap_image mirror(const bitmap_image &source, char direction) {
    bitmap_image res(source.width(), source.height());
    for (unsigned int i = 0; i < source.width(); ++i) {
        for (unsigned int j = 0; j < source.height(); ++j) {
            switch (direction) {
                case 'h':
                    res.set_pixel(source.width() - i - 1, j, source.get_pixel(i, j));
                    break;
                case 'w':
                    res.set_pixel(i, source.height() - j - 1, source.get_pixel(i, j));
                    break;
                default:
                    std::cout << "default case reached! Something went wrong..." << std::endl;
                    break;
            }
        }
    }
    return res;
}

bitmap_image rotate_ccw(const bitmap_image &source) {
    bitmap_image res(source.height(), source.width());

    for (unsigned int i = 0; i < source.width(); ++i) {
        for (unsigned int j = 0; j < source.height(); ++j) {
            res.set_pixel(j, source.width() - 1 - i, source.get_pixel(i, j));
        }
    }

    return res;
}

bitmap_image leave_single_color(const bitmap_image &source, char color) {
    bitmap_image res(source.width(), source.height());
    rgb_t pixel = make_colour(0, 0, 0);
    uint8_t yuv_state;
    for (unsigned int i = 0; i < source.width(); ++i) {
        for (unsigned int j = 0; j < source.height(); ++j) {
            //switch
            if (color == 'r') {
                pixel = source.get_pixel(i, j);
                pixel.blue = 0;
                pixel.green = 0;
            } else if (color == 'g') {
                pixel = source.get_pixel(i, j);
                pixel.blue = 0;
                pixel.red = 0;
            } else if (color == 'b') {
                pixel = source.get_pixel(i, j);
                pixel.red = 0;
                pixel.green = 0;
            } else if (color == 'y') {
                pixel = source.get_pixel(i, j);
                yuv_state = static_cast<uint8_t>(0.257 * pixel.red + 0.504 * pixel.green + 0.098 * pixel.blue + 16);
                pixel.blue = yuv_state;
                pixel.green = yuv_state;
                pixel.red = yuv_state;
            } else if (color == 'B') {
                pixel = source.get_pixel(i, j);
                yuv_state = static_cast<uint8_t>(-0.148 * pixel.red - 0.291 * pixel.green + 0.439 * pixel.blue + 128);
                pixel.blue = yuv_state;
                pixel.green = yuv_state;
                pixel.red = yuv_state;
            } else if (color == 'R') {
                pixel = source.get_pixel(i, j);
                yuv_state = static_cast<uint8_t>(0.439 * pixel.red - 0.368 * pixel.green - 0.071 * pixel.blue + 128);
                pixel.blue = yuv_state;
                pixel.green = yuv_state;
                pixel.red = yuv_state;
            }
            res.set_pixel(i, j, pixel);
        }
    }
    return res;
}

bitmap_image generate_rgb_back(const bitmap_image &y, const bitmap_image &cb, const bitmap_image &cr) {
    bitmap_image res(y.width(), y.height());
    rgb_t pixel;
    for (unsigned int i = 0; i < res.width(); ++i) {
        for (unsigned int j = 0; j < res.height(); ++j) {
            pixel.red = static_cast<unsigned char>(1.164 * (y.get_pixel(i, j).green - 16) +
                                                   1.596 * (cr.get_pixel(i, j).green - 128));
            pixel.green = static_cast<unsigned char>(1.164 * (y.get_pixel(i, j).green - 16) -
                                                     0.813 * (cr.get_pixel(i, j).green - 128) -
                                                     0.392 * (cb.get_pixel(i, j).green - 128));
            pixel.blue = static_cast<unsigned char>(1.164 * (y.get_pixel(i, j).green - 16) +
                                                    2.017 * (cb.get_pixel(i, j).green - 128));
            res.set_pixel(i, j, pixel);
        }
    }
    return res;
}

//14
std::vector<std::vector<std::vector<int>>> DPCM(const bitmap_image &source, char channel) {

    std::vector<std::vector<std::vector<int>>> res(4, std::vector<std::vector<int>>(source.width(),
                                                                                    std::vector<int>(source.height())));
    for (int k = 0; k < 4; ++k) {
        res[k].reserve(source.width() * source.height());
        for (int i = 1; i < source.width(); i++) {
            res[k][i].reserve(source.height());
        }
    }
    for (int i = 1; i < source.width(); i++) {
        for (int j = 1; j < source.height(); j++) {
            switch (channel) {
                case 'r':
                    res[0][i - 1][j - 1] = source.get_pixel(i, j).red - source.get_pixel(i, j - 1).red + 128;
                    res[1][i - 1][j - 1] = source.get_pixel(i, j).red - source.get_pixel(i - 1, j).red + 128;
                    res[2][i - 1][j - 1] = source.get_pixel(i, j).red - source.get_pixel(i - 1, j - 1).red + 128;
                    res[3][i - 1][j - 1] = source.get_pixel(i, j).red -
                                           (source.get_pixel(i - 1, j - 1).red + source.get_pixel(i, j - 1).red +
                                            source.get_pixel(i - 1, j).red) / 3 + 128;
                    break;
                case 'g':
                    res[0][i - 1][j - 1] = source.get_pixel(i, j).green - source.get_pixel(i, j - 1).green + 128;
                    res[1][i - 1][j - 1] = source.get_pixel(i, j).green - source.get_pixel(i - 1, j).green + 128;
                    res[2][i - 1][j - 1] = source.get_pixel(i, j).green - source.get_pixel(i - 1, j - 1).green + 128;
                    res[3][i - 1][j - 1] = source.get_pixel(i, j).green -
                                           (source.get_pixel(i - 1, j - 1).green + source.get_pixel(i, j - 1).green +
                                            source.get_pixel(i - 1, j).green) / 3 + 128;
                    break;
                case 'b':
                    res[0][i - 1][j - 1] = source.get_pixel(i, j).blue - source.get_pixel(i, j - 1).blue + 128;
                    res[1][i - 1][j - 1] = source.get_pixel(i, j).blue - source.get_pixel(i - 1, j).blue + 128;
                    res[2][i - 1][j - 1] = source.get_pixel(i, j).blue - source.get_pixel(i - 1, j - 1).blue + 128;
                    res[3][i - 1][j - 1] = source.get_pixel(i, j).blue -
                                           (source.get_pixel(i - 1, j - 1).blue + source.get_pixel(i, j - 1).blue +
                                            source.get_pixel(i - 1, j).blue) / 3 + 128;
                    break;
                default:
                    std::cout << "default case reached! Something went wrong..." << std::endl;
                    break;
            }
        }
    }
    return res;
}

//16
std::vector<double>
dpcm_enthropies(const std::vector<std::vector<std::vector<int>>> &data, const bitmap_image &source) {
    std::vector<double> res(256);
    std::vector<double> results(4);
    for (int k = 0; k < 4; k++) {
        double result = 0;
        for (int i = 0; i < source.width() - 1; i++) {
            for (int j = 0; j < source.height() - 1; j++) {
                res[abs(data[k][i][j]) % 256]++;
            }
        }

        for (int j = 0; j < 256; j++) {
            res[j] /= source.height() * source.width();
            if (res[j] != 0) {
                result -= res[j] * log2(res[j]);
            }
            res[j] = 0;
        }
        results[k] = result;
    }
    return results;
}