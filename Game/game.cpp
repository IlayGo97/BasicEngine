#include "game.h"
#include "stb_image.h"
#include <iostream>
#include <glm/gtc/matrix_transform.hpp>
#include <numbers>


static void printMat(const glm::mat4 mat) {
    std::cout << " matrix:" << std::endl;
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++)
            std::cout << mat[j][i] << " ";
        std::cout << std::endl;
    }
}

Game::Game() : Scene() {
}

Game::Game(float angle, float relationWH, float near1, float far1) : Scene(angle, relationWH, near1, far1) {
}

#define color_size_bytes 4
#define black(ARR, X, Y) ARR[X][(Y)] = 0; ARR[X][(Y)+1] = 0; ARR[X][(Y)+2] = 0; ARR[X][(Y)+3] = 255;
#define white(ARR, X, Y) ARR[X][(Y)] = 255; ARR[X][(Y)+1] = 255; ARR[X][(Y)+2] = 255; ARR[X][(Y)+3] = 255;
#define set_val(ARR, X, Y, val) ARR[X][(Y)] = val; ARR[X][(Y)+1] = val; ARR[X][(Y)+2] = val; ARR[X][(Y)+3] = 255;
#define set_val(ARR, INDEX, val) ARR[INDEX] = val; ARR[INDEX+1] = val; ARR[INDEX+2] = val; ARR[INDEX+3] = 255;
#define to_index(i, j) i * width * color_size_bytes + j * color_size_bytes
#define to_index_normal(i, j) (i) * width + (j)
#define pixel_average(ARR, X, Y) ((ARR[X][(Y)] + ARR[X][(Y)+1] + ARR[X][(Y)+2])/3)
#define PI 3.141592654

static unsigned char **single_array_to_multi(unsigned char *data, int width, int height) {
    unsigned char **output = (unsigned char **) malloc(height * sizeof(unsigned char *));
    for (int i = 0; i < height; i++)
        output[i] = (unsigned char *) malloc(width * 4);
    int data_size = color_size_bytes * width * height;
    for (int i = 0; i < data_size; i++) {
        output[i / (width * 4)][i % (4 * width)] = data[i];
    }
    return output;
}

static unsigned char *TwoD2OneD(unsigned char **arr, int width, int height) {
    unsigned char *output = (unsigned char *) malloc(width * height * sizeof(uint64_t));
    for (int i = 0; i < height; i++)
        for (int j = 0; j < width * 4; j++)
            output[i * width * 4 + j] = arr[i][j];
    return output;
}

static double convert_to_degree(int rad) {
    return rad * (180 / PI) + 180;
}

static int round(int val, int val_amount) {
    return (val / val_amount) * val_amount;
}

static void inc(unsigned char *data, unsigned char amount, int i, int j, int width, int height) {
    int val = data[to_index(i, j)];
    val = val + amount;
    set_val(data, to_index(i, j), val);
}

void Game::AddFloydSteinbergText() {
    int width, height, numComponents;
    double alpha = (double) 7 / 16;
    double beta = (double) 3 / 16;
    double gamma = (double) 5 / 16;
    double delta = (double) 1 / 16;
    unsigned char *data = stbi_load("../res/textures/lena256.jpg", &width, &height, &numComponents, 4);
    auto *output = (unsigned char *) malloc(height * width * color_size_bytes * sizeof(unsigned char));
    memcpy(output, data, height * width * color_size_bytes * sizeof(unsigned char));
    for (int i = 0; i < height - 1; i++)
        for (int j = 0; j < width - 1; j++) {
            int P = round(data[to_index(i, j)], 16);
            int e = data[to_index(i, j)] - P;
            set_val(output, to_index(i, j), P)
            inc(data, alpha * e, i, j + 1, width, height);
            inc(data, beta * e, i + 1, j - 1, width, height);
            inc(data, gamma * e, i + 1, j, width, height);
            inc(data, delta * e, i + 1, j + 1, width, height);
        }
    AddTexture(width, height, output);
}

void Game::AddHalftonePatternText() {
    int width, height, numComponents;
    unsigned char *data = stbi_load("../res/textures/lena256.jpg", &width, &height, &numComponents, 4);
    unsigned char **data2d = single_array_to_multi(data, width, height);
    unsigned char **output2d = (unsigned char **) malloc(2 * height * sizeof(unsigned char *));
    for (int i = 0; i < height * 2; i++)
        output2d[i] = (unsigned char *) malloc(2 * width * color_size_bytes * sizeof(unsigned char));
    for (int i = 0; i < height; i++)
        for (int j = 0; j < width * color_size_bytes; j = j + color_size_bytes) {
            int average = pixel_average(data2d, i, j);
            if (average < 0.20 * 255) {
                black(output2d, 2 * i, 2 * j)
                black(output2d, 2 * i, 2 * j + color_size_bytes)
                black(output2d, 2 * i + 1, 2 * j)
                black(output2d, 2 * i + 1, 2 * j + color_size_bytes)
            } else if (average < 0.4 * 255) {
                white(output2d, 2 * i, 2 * j)
                black(output2d, 2 * i, 2 * j + color_size_bytes)
                black(output2d, 2 * i + 1, 2 * j)
                black(output2d, 2 * i + 1, 2 * j + color_size_bytes)
            } else if (average < 0.6 * 255) {
                white(output2d, 2 * i, 2 * j)
                black(output2d, 2 * i, 2 * j + color_size_bytes)
                black(output2d, 2 * i + 1, 2 * j)
                white(output2d, 2 * i + 1, 2 * j + color_size_bytes)
            } else if (average < 0.8 * 255) {
                white(output2d, 2 * i, 2 * j)
                white(output2d, 2 * i, 2 * j + color_size_bytes)
                black(output2d, 2 * i + 1, 2 * j)
                white(output2d, 2 * i + 1, 2 * j + color_size_bytes)
            } else {
                white(output2d, 2 * i, 2 * j);
                white(output2d, 2 * i, 2 * j + color_size_bytes);
                white(output2d, 2 * i + 1, 2 * j);
                white(output2d, 2 * i + 1, 2 * j + color_size_bytes);
            }
        }
    AddTexture(2 * width, 2 * height, TwoD2OneD(output2d, width * 2, height * 2));
}

//void Game::AddHalftonePatternText() {
//    int width, height, numComponents;
//    unsigned char* data = stbi_load("../res/textures/lena256.jpg", &width, &height, &numComponents, 4);
//    unsigned char** data2d = single_array_to_multi(data, width, height);
//    unsigned char** output2d = (unsigned  char**)malloc(height * sizeof (unsigned  char*));
//    for(int i = 0; i < height; i ++)
//        output2d[i] = (unsigned  char*)malloc( width * color_size_bytes * sizeof (unsigned char));
//    for(int i = 0; i < height - 1; i = i + 2)
//        for(int j = 0; j < (width - 1) * color_size_bytes; j = j + 2 * color_size_bytes)
//        {
//            int average = pixel_average(data2d,i,j) + pixel_average(data2d,i,j + 1) + pixel_average(data2d,i + 1,j) + pixel_average(data2d,i + 1,j + 1);
//            average = average /  4;
//            if (average < 0.20 * 255){
//                black(output2d,i,j)
//                black(output2d,i,j + color_size_bytes)
//                black(output2d,i + 1,j)
//                black(output2d,i + 1,j + color_size_bytes)
//            } else if (average < 0.4 * 255) {
//                white(output2d,i,j)
//                black(output2d,i,j + color_size_bytes)
//                black(output2d,i + 1,j)
//                black(output2d,i + 1,j + color_size_bytes)
//            } else if (average < 0.6 * 255) {
//                white(output2d,i,j)
//                black(output2d,i,j + color_size_bytes)
//                black(output2d,i + 1,j)
//                white(output2d,i + 1,j + color_size_bytes)
//            } else if (average < 0.8 * 255) {
//                white(output2d,i,j)
//                white(output2d,i,j + color_size_bytes)
//                black(output2d,i + 1,j)
//                white(output2d,i + 1,j + color_size_bytes)
//            } else{
//                white(output2d,i,j)
//                white(output2d,i,j + color_size_bytes)
//                white(output2d,i + 1,j)
//                white(output2d,i + 1,j + color_size_bytes)
//            }
//        }
//    AddTexture( width,height, TwoD2OneD(output2d, width, height));
//}



unsigned char *
applyConv(const unsigned char *data, int width, int height, const glm::detail::tmat3x3<float, glm::highp> &filter);

int *applyGrad(const unsigned char *data, int width, int height, const glm::detail::tmat3x3<float, glm::highp> &filter);

unsigned char *applyGaus3(unsigned char *data, int width, int height) {
    auto filter = glm::detail::tmat3x3<float, glm::highp>(
            1, 2, 1,
            2, 4, 2,
            1, 2, 1);
    float sum = 16.0;
    filter = filter / sum;
    unsigned char *output = applyConv(data, width, height, filter);
    return output;
}

int *applyGradx(unsigned char *data, int width, int height) {
    auto filter = glm::detail::tmat3x3<float, glm::highp>(
            -1, 0, 1,
            -2, 0, 2,
            -1, 0, 1);
    int *output = applyGrad(data, width, height, filter);
    return output;
}

int *
applyGrad(const unsigned char *data, int width, int height, const glm::detail::tmat3x3<float, glm::highp> &filter) {
    int *output = (int *) malloc(width * height * sizeof(int));
    for (int i = 1; i < height - 1; i++)
        for (int j = 1; j < width - 1; j++) {
            output[to_index_normal(i, j)] = 0;
            for (int h = i; h < i + 3; h++)
                for (int w = j; w < j + 3; w++) {
                    int newval = output[to_index_normal(i, j)] + filter[h - i][w - j] * data[to_index(h, w)];
                    output[to_index_normal(i, j)] = newval;
                }
        }
    return output;
}

int *applyGrady(unsigned char *data, int width, int height) {
    auto filter = glm::detail::tmat3x3<float, glm::highp>(
            -1, -2, -1,
            0, 0, 0,
            1, 2, 1);
    int *output = applyGrad(data, width, height, filter);
    return output;
}

unsigned char *
applyConv(const unsigned char *data, int width, int height, const glm::detail::tmat3x3<float, glm::highp> &filter) {
    auto *output = (unsigned char *) malloc(height * width * color_size_bytes * sizeof(unsigned char));
    memcpy(output, data, height * width * color_size_bytes * sizeof(unsigned char));
    for (int i = 1; i < height - 1; i++)
        for (int j = 1; j < width - 1; j++) {
            set_val(output, to_index(i, j), 0);
            for (int h = i; h < i + 3; h++)
                for (int w = j; w < j + 3; w++) {
                    int newval = output[to_index(i, j)] + filter[h - i][w - j] * data[to_index(h, w)];
                    newval = std::max(0, newval);
                    newval = std::min(255, newval);
                    set_val(output, to_index(i, j), newval);
                }
        }
    return output;
}

int *
getGrad(int *gradX, int *gradY, int height, int width) {
    auto *output = (int *) malloc(height * width * sizeof(int));
//    int count = 0;
    memset(output, 0, height * width * sizeof(int));
    for (int i = 1; i < height - 1; i++)
        for (int j = 1; j < width - 1; j++) {
            int index = to_index_normal(i, j);
            int new_val = gradY[index] * gradY[index] + gradX[index] * gradX[index];
            new_val = std::sqrt(new_val);
            output[to_index_normal(i, j)] = new_val;
//            count++;
        }
//    std::cout << count << std::endl;
    return output;
}

unsigned char *GradientOrientation(int *grad_x, int *grad_y, int height, int width) {
    auto *gradient_orientation = (unsigned char *) malloc (height * width * sizeof(unsigned char));
    for (int y = 1; y < height - 1; y++) {
        for (int x = 1; x < width - 1; x++) {
            int index = to_index_normal(y, x);
            unsigned char theta = convert_to_degree(atan2(grad_y[index], grad_x[index]));
            gradient_orientation[index] = theta;
        }
    }
}

unsigned char *NonMaximumSuppression(int *gradient_magnitude, unsigned char *gradient_direction, int height, int width) {
    auto *output = (unsigned char *) calloc (height * width * 4 * sizeof(unsigned char), sizeof(unsigned char));
    for (int y = 1; y < height - 1; y++) {
        for (int x = 1; x < width - 1; x++) {
            int index = to_index_normal(y, x);
            int before_pixel = 0;
            int after_pixel = 0;
            unsigned char direction = gradient_direction[index];
            std::cout << direction << std::endl;
            if ((0 <= direction < PI / 8) || (15 * PI / 8 <= direction <= 2 * PI)) {
                before_pixel = gradient_magnitude[to_index_normal(y, x - 1)];
                after_pixel = gradient_magnitude[to_index_normal(y, x + 1)];
            } else if ((PI / 8 <= direction * 3 * PI / 8) || (9 * PI / 8 <= direction < 11 * PI / 8)) {
                before_pixel = gradient_magnitude[to_index_normal(y + 1, x - 1)];
                after_pixel = gradient_magnitude[to_index_normal(y - 1, x + 1)];
            } else if ((3 * PI / 8 <= direction < 5 * PI / 8) || (11 * PI / 8 <= direction < 13 * PI / 8)) {
                before_pixel = gradient_magnitude[to_index_normal(y - 1, x)];
                after_pixel = gradient_magnitude[to_index_normal(y + 1, x)];
            } else {
                before_pixel = gradient_magnitude[to_index_normal(y - 1, x - 1)];
                after_pixel = gradient_magnitude[to_index_normal(y + 1, x + 1)];

            }
            if (gradient_magnitude[to_index_normal(y, x)] >= before_pixel &&
                gradient_magnitude[to_index_normal(y, x)] >= after_pixel) {
                output[to_index(y, x)] = gradient_magnitude[to_index_normal(y, x)];
            }
        }
    }
    return output;
}


unsigned char *
Thresholding(int *grad, int highthresh, int lowthresh, int height, int width) {
    auto *output = (unsigned char *) malloc(height * width * color_size_bytes * sizeof(unsigned char));
    memset(output, 0, height * width * color_size_bytes * sizeof(unsigned char));
    for (int i = 1; i < height - 1; i++)
        for (int j = 1; j < width - 1; j++) {
            int index = to_index(i, j);
            int val = grad[to_index_normal(i, j)];
            if (val < lowthresh) {
                set_val(output, index, 0)
            } else if (val < highthresh) {
                set_val(output, index, 0)
            } else {
                set_val(output, index, 255)
            }
        }
    return output;
}

void Game::AddEdgesText() {
    int max_grad = std::sqrt(255 * 255 * 2);
    int width, height, numComponents;
    unsigned char *data = stbi_load("../res/textures/lena256.jpg", &width, &height, &numComponents, 4);
    data = applyGaus3(data, width, height);
    int *gradX = applyGradx(data, width, height);
    int *gradY = applyGrady(data, width, height);
    int *grad = getGrad(gradX, gradY, height, width);
    unsigned char *theta = GradientOrientation(gradX, gradY, height, width);
    data = NonMaximumSuppression(grad, theta, height, width);
    data = Thresholding(grad, max_grad * 0.5, max_grad * 0.3, height, width);
    delete[]  gradX;
    delete[] gradY;
    delete[] grad;
    AddTexture(width, height, data);
}

void Game::Init() {

    AddShader("../res/shaders/pickingShader");
    AddShader("../res/shaders/basicShader");
    AddTexture("../res/textures/lena256.jpg", false); // 0
    AddEdgesText();                                                     // 1
    AddHalftonePatternText();                                           // 2
    AddFloydSteinbergText();                                            // 3
    AddShape(Plane, -1, TRIANGLES);
    AddShape(Plane, -1, TRIANGLES);
    AddShape(Plane, -1, TRIANGLES);
    AddShape(Plane, -1, TRIANGLES);
    for (int i = 0; i < 4; i++) {
        pickedShape = i;
        SetShapeTex(i, i);
    }
    MoveCamera(0, zTranslate, 10);
    pickedShape = -1;

    //ReadPixel(); //uncomment when you are reading from the z-buffer
}

void Game::Update(const glm::mat4 &MVP, const glm::mat4 &Model, const int shaderIndx) {
    Shader *s = shaders[shaderIndx];
    int r = ((pickedShape + 1) & 0x000000FF) >> 0;
    int g = ((pickedShape + 1) & 0x0000FF00) >> 8;
    int b = ((pickedShape + 1) & 0x00FF0000) >> 16;
    s->Bind();
    s->SetUniformMat4f("MVP", MVP);
    s->SetUniformMat4f("Normal", Model);
    s->SetUniform4f("lightDirection", 0.0f, 0.0f, -1.0f, 0.0f);
    if (shaderIndx == 0)
        s->SetUniform4f("lightColor", r / 255.0f, g / 255.0f, b / 255.0f, 1.0f);
    else
        s->SetUniform4f("lightColor", 0.7f, 0.8f, 0.1f, 1.0f);
    s->Unbind();
}

void Game::WhenRotate() {
}

void Game::WhenTranslate() {
}

void Game::Motion() {
    if (isActive) {
    }
}

Game::~Game(void) {
}
