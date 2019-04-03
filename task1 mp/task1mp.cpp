#include <iostream>
#include<vector>
#include <iostream>
#include <vector>
#include <fstream>
#include <numeric>
#include <chrono>
#include <fstream> // To be removed later
#include<omp.h>
#include<stdio.h>


// Using lodepng as png image manipulator
#include "lodepng.h" // PNG loader 
//#define N 4  // Number of threads

using namespace std;
/*
	---------------------------------------- Global variables -------------------------------------------------
*/
const char *FILE_1 = "Images/im0.png"; // First image file name
const char *FILE_2 = "Images/im1.png"; // Second image file name
const char *FILE_D1 = "Images/d_im0.png"; // Downsampled first image file name
const char *FILE_D2 = "Images/d_im1.png"; // Downsampled second image file name
const char *FILE_L2R = "Images/left_on_right.png";
const char *FILE_R2L = "Images/right_on_left.png";
const char *FILE_CROSS_CHECKED = "Images/cross_checked.png";
const char *FILE_OCCLUSION_FILLED = "Images/occlusion_filled.png";
const int M = 16; // Downsampling factor
int WIN_SIZE = 9;//windowsize
int N = 4;
unsigned MAX_DISP = (unsigned)260 / 4; // Initial value provided in "calib.txt"
const int THRESHOLD = 8 * 256 / 65; // Threshold used on cross check (mapped)
vector<std::vector<unsigned char> > final_img; // Final output image. Assigned only at the end of the program

struct img_manipulator
{
	vector<unsigned char> image; //the raw pixels	// Freed once gray2D ready 
	vector<std::vector<unsigned char> > gray2D;
	vector<std::vector<unsigned char> > disp_img; // It will contain disp image when depthEstim is called

	unsigned width, height;
};

/*
	---------------------------------------- Function signature -------------------------------------------------
*/

void loadImages(vector<unsigned char>& image, unsigned& width, unsigned& height, const char* filename);
void writeImage(vector<std::vector<unsigned char> >& image, const char* filename, bool color);
void rgba2gray(vector<unsigned char>& image, vector<vector<unsigned char>>& gray2D, unsigned& width, unsigned& height, int M);











/*
double ZNCC_formula(unsigned x, unsigned y, unsigned b, vector<std::vector<unsigned char> >ilf, vector<std::vector<unsigned char> > irf, float war, float wal, int d) {

	
	std::vector<std::vector<double>> il;
	il.reserve(ilf.size());
	for (auto&& v : ilf) il.emplace_back(std::begin(v), std::end(v));
	std::vector<std::vector<double>> ir;
	ir.reserve(irf.size());
	for (auto&& v : irf) ir.emplace_back(std::begin(v), std::end(v));

	double k = 0, p = 0, m = 0 ,ZNCC=0;
	
	for (int j = y ; j <b + y-1 ; j++) {
		for (int i = x + d; i < b + x + d-1; i++) {
			k = k + (il[j][i] - wal)*(ir[j][i-d] - war);
			p = p + pow((il[j][i] - wal), 2);
			m = m + pow((ir[j][i - d] - war), 2);
		}


	}
	/*

	for (int j = y ; j < b + y-1 ; j++) {
		for (int i = x + d; i < b + x + d-1; i++) {
			p = p+ pow((il[j][i]  - wal),2);

		}


	}


	for (int j = y ; j < b + y-1 ; j++) {
		for (int i = x + d; i< b + x + d-1; i++) {
			m = m + pow((ir[j][i-d] - war),2);

		}


	}
	
	ZNCC = k / (sqrt(p)*sqrt(m));
	
	
	return ZNCC; 



}

*/




int main()
{
	


	// // Loading images to memory
	// Creating an image manipulator structure object
	img_manipulator img1;
	img_manipulator img2;

	// Let's load the images into memory
	loadImages(img1.image, img1.width, img1.height, FILE_1); // Loading first image
	loadImages(img2.image, img2.width, img2.height, FILE_2); // Loading second image

	// Print out inital loaded image characteristics
	std:: cout << "Height of img1: " << img1.height << " Width of img1: " << img1.width << endl;
	std:: cout << "Image: " << img1.image.size() << endl;

	// Gray2D
	std::cout << "Applying gray on image1... " << endl;
	rgba2gray(img1.image, img1.gray2D, img1.width, img1.height, M);
	std::cout << "Applying gray on image2... " << endl;
	rgba2gray(img2.image, img2.gray2D, img2.width, img2.height, M);

	// Print out inital loaded image characteristics
	std::cout << "Height of gray2D 1 " << img1.gray2D.size() << " Width of img1: " << img1.gray2D[0].size() << endl;
	std::cout << "Height of gray2D 2  " << img2.gray2D.size() << " Width of img2: " << img2.gray2D[0].size() << endl;
	
	double sum_1 = 0, sum_2 = 0, sum_3 = 0, sum_4 = 0, window_mean_1=0, window_mean_2=0,  current_max_sum=-2, current_max_sum__2=-2, window_mean_3 = 0, window_mean_4 = 0, current_max_sum_2 = -2;
	int exr=0;
	unsigned best_disparity=0, best_disparity_2 = 0;
	double k = 0, p = 0, m = 0, Zncc_value = -1, k_1= 0, p_1 = 0, m_1 = 0, Zncc_value_2 = -1;
	vector<std::vector<unsigned char> > vect;
	vector<std::vector<unsigned char> > vect_2;
	
	
	std::vector<std::vector<double>> il;
	il.reserve(img1.gray2D.size());
	for (auto&& v : img1.gray2D) il.emplace_back(std::begin(v), std::end(v));
	std::vector<std::vector<double>> ir;
	ir.reserve(img2.gray2D.size());
	for (auto&& v : img2.gray2D) ir.emplace_back(std::begin(v), std::end(v));
	
//ZNCC Algorithm img1.height-9	img1.width-9
#pragma omp parallel for default(none)
	for (unsigned j = 0; j < img1.height - WIN_SIZE; j++) {
		vector<unsigned char> temp;
		
		//cout << img1.height << endl;
		for (unsigned i = 0; i < img1.width -WIN_SIZE; i++) {
			current_max_sum = -2;
			

			//best_disparity = 0;
			//cout <<i << endl;

			for (unsigned d = 0; d < 64; d++) {
				sum_1 = 0;
				sum_2 = 0;
				
				k = 0;
				
				p = 0; 
				m = 0;
				Zncc_value=0;
				Zncc_value_2 = 0;
				for (unsigned Win_Y = j; Win_Y < WIN_SIZE +j; Win_Y++) {

					for (unsigned Win_X = i+d; Win_X < WIN_SIZE +i+d; Win_X++) {
						if ((Win_Y + 8 >= img1.height - WIN_SIZE) || ((Win_X + 8 >= img1.width - WIN_SIZE))) {
							exr = 1;
							break;
						}
						
						//calculate the Mean Value for each window 
						
						sum_1 += img1.gray2D[Win_Y][Win_X];
						sum_2 += img2.gray2D[Win_Y][Win_X-d];
						
					}
					if (exr == 1) {
						
						break;
					}
				}
				window_mean_1 = sum_1 /( WIN_SIZE* WIN_SIZE);
				window_mean_2 = sum_2 / (WIN_SIZE* WIN_SIZE);
				
				
				
				
				if (exr==1) {
					exr = 0;
					break;
					}
						//calcualte the ZNCC Value for each Window
				//Zncc_value = ZNCC_formula(i, j, 9, img1.gray2D, img2.gray2D, window_mean_1, window_mean_2, d);


				
				for (int zj = j; zj < WIN_SIZE + j; zj++) {
					
					 for (int zi = i + d; zi < WIN_SIZE + i + d; zi++) {
						 k += (il[zj][zi] - window_mean_1)*(ir[zj][zi - d] - window_mean_2);
						 p += pow((il[zj][zi] - window_mean_1), 2);
						 m += pow((ir[zj][zi - d] - window_mean_2), 2);

						 
					 }
				 }
				
				 Zncc_value = k / (sqrt(p)*sqrt(m));
				

				
				//cout <<d << endl;
				//cout << Zncc_value << endl;
				if (Zncc_value > current_max_sum) {

					//update current_max_sum
					current_max_sum = Zncc_value;

					//update Best Disparity Value
					best_disparity =( d*255)/64;
				
				}



				
				//cout << d << endl;
			
				//cout << current_max_sum << endl;
			}

		//Disparity_image_pixel=Best_Disparity_Value
			//cout << best_disparity << endl;
			temp.push_back( best_disparity);
			
			//cout << best_disparity << endl;
			
			
		
		}
		

		vect.push_back(temp);
		

		//std:: cout << best_disparity << endl;
		//cout << j << endl;
	}
	//zncc2
	exr = 0;
#pragma omp parallel for collapse(2)  num_threads(4) schedule(static, 3)
	for (unsigned j = 0; j < img1.height - WIN_SIZE; j++) {
		
		vector<unsigned char> temp_2;
		//cout << img1.height << endl;
		for (unsigned i = 0; i < img1.width - WIN_SIZE; i++) {
			
			current_max_sum_2 = -2;

			//best_disparity = 0;
			//cout <<i << endl;
			for (unsigned d = 0; d < 64; d++) {
				
				sum_3 = 0;
				sum_4 = 0;
				
				k_1 = 0;
				p_1 = 0;
				m_1 = 0;
				
				
				for (unsigned Win_Y = j; Win_Y < WIN_SIZE + j; Win_Y++) {

					for (unsigned Win_X = i - d; Win_X < WIN_SIZE + i - d; Win_X++) {
						if ((Win_Y + WIN_SIZE-1 >= img1.height - WIN_SIZE) || ((Win_X + WIN_SIZE-1 >= img1.width - WIN_SIZE))) {
							exr = 1;
							break;
						}
						if ((Win_X  < 0) ){
							exr = 1;
							break;
						}
						//calculate the Mean Value for each window 

						
						sum_3 += img1.gray2D[Win_Y][Win_X+d];
						sum_4 += img2.gray2D[Win_Y][Win_X ];
					}
					if (exr == 1) {

						break;
					}
				}
				
				window_mean_3 = sum_3 / (WIN_SIZE* WIN_SIZE);
				window_mean_4 = sum_4 / (WIN_SIZE* WIN_SIZE);



				if (exr == 1) {
					exr = 0;
					break;
				}
				//calcualte the ZNCC Value for each Window
		//Zncc_value = ZNCC_formula(i, j, 9, img1.gray2D, img2.gray2D, window_mean_1, window_mean_2, d);

#pragma omp parallel for collapse(2)   num_threads(4) schedule(static, 3)
				for (int kj = j; kj < WIN_SIZE + j; kj++) {

					for (int ki = i - d; ki < WIN_SIZE + i - d; ki++) {


						k_1 += (ir[kj][ki] - window_mean_4)*(il[kj][ki + d] - window_mean_3);
						p_1 += pow((ir[kj][ki] - window_mean_4), 2);
						m_1 += pow((il[kj][ki + d] - window_mean_3), 2);
					}
				}
				
				Zncc_value_2 = k_1 / (sqrt(p_1)*sqrt(m_1));





				if (Zncc_value_2 > current_max_sum_2) {

					//update current_max_sum
					current_max_sum_2 = Zncc_value_2;

					//update Best Disparity Value
					best_disparity_2 = (d * 255) / 64;

				}
				//cout << d << endl;

				//cout << current_max_sum << endl;
			}

			//Disparity_image_pixel=Best_Disparity_Value
				//cout << best_disparity << endl;
			
			temp_2.push_back(best_disparity_2);
			//cout << best_disparity << endl;



		}


		
		vect_2.push_back(temp_2);

		//std:: cout << best_disparity << endl;
		//cout << j << endl;
	}
	std::vector<std::vector<unsigned char>> img_disp;
	img_disp.reserve(vect.size());
	for (auto&& v : vect) img_disp.emplace_back(std::begin(v), std::end(v));
	writeImage(img_disp, FILE_L2R, false);

	std::vector<std::vector<unsigned char>> img_disp_2;
	img_disp_2.reserve(vect_2.size());
	for (auto&& v : vect_2) img_disp_2.emplace_back(std::begin(v), std::end(v));
	writeImage(img_disp_2, FILE_R2L, false);

	//Cross Checking 

	vector<std::vector<unsigned char> > vect_3;
	const int th = 25;
	double diff;
#pragma omp parallel for collapse(2)   schedule(static, 3)
	for (unsigned j = 0; j < img1.height- WIN_SIZE; j++) {

		vector<unsigned char> temp_3;
	
		for (unsigned i = 0; i < img1.width- WIN_SIZE; i++) {
			diff = abs(vect[j][i] - vect_2[j][i]);
			if (diff >= th) {
				temp_3.push_back(0);
			}
			else {
				temp_3.push_back((( vect_2[j][i])));
			}

		}

		vect_3.push_back(temp_3);
	}


	std::vector<std::vector<unsigned char>> img_disp_3;
	img_disp_3.reserve(vect_3.size());
	for (auto&& v : vect_3) img_disp_3.emplace_back(std::begin(v), std::end(v));
	writeImage(img_disp_3, FILE_CROSS_CHECKED, false);
	

	//occlusion filling 
	
#pragma omp parallel for collapse(2)   schedule(static, 3)
	for (unsigned j = 0; j < img1.height- WIN_SIZE; j++) {

		vector<unsigned char> temp_3;

		for (unsigned i = 0; i < img1.width- WIN_SIZE; i++) {
			
			if (vect_3[j][i] == 0) {
				for (unsigned pj = 1; pj < 200; pj++) {
					
					double a[8] = { 0,0,0,0,0,0,0,0 }, b = 0, c = 0;
					if ((i + pj) < (img1.width-9)) {
						a[3] = vect_3[j][i + pj];
					}
					
					/*
					if ((i - pj) > 0) {
						a[7] = vect_3[j][i - pj];
					}
					*/
					if ((j + pj)  < (img1.height-9)) {
						a[1] = vect_3[j + pj][i];
					}
					//cout << j << endl;
					/*
					if ((j - pj) > 0) {
						
						a[5] = vect_3[j - pj][i];
					}
					*/
					if (((j + pj) < (img1.height-9))&((i + pj) < (img1.width-9))) {
						a[2] = vect_3[j + pj][i + pj];
					}
					/*
					if (((i + pj) < (img1.width-9))&((i - pj) > 0)) {
						a[0] = vect_3[j + pj][i - pj];
					}
					*/
					/*
					if (((j - pj) > 0)&((i + pj) < (img1.width-9))) {
						a[4] = vect_3[j - pj][i + pj];
					}
					*/
					/*
					if (((j - pj) > 0)&((i - pj) > 0)) {
						a[6] = vect_3[j - pj][i - pj];
					}
					*/
					
					double tempy = 0;
					for (unsigned oj = 0; oj < 8; oj++) {
						if (a[oj] > 0) {
							b = b + 1;

						}
					}
					
					
					for (unsigned aj = 0; aj < 8; aj++) {
						
						
						
						
						if (a[aj] > tempy) {
							//b = b + 1;
							tempy = a[aj];
						}

					}
					if (b >= 2) {
						vect_3[j][i] = tempy;
						break;
					}

				}




			}
			
			}

		}
	std::vector<std::vector<unsigned char>> img_disp_4;
	img_disp_4.reserve(vect_3.size());
	for (auto&& v : vect_3) img_disp_4.emplace_back(std::begin(v), std::end(v));
	writeImage(img_disp_4, FILE_OCCLUSION_FILLED, false);
		
	
	
	
	
	system("pause");
	return 0;
}


// Source: https://raw.githubusercontent.com/lvandeve/lodepng/master/examples/example_decode.cpp
//Load PNG file from disk to memory first, then decode to raw pixels in memory.
void loadImages(vector<unsigned char>& image, unsigned& width, unsigned& height, const char* filename) {
	std::vector<unsigned char> png; // Image will be represented in 8-bit unsinged representation

	//load and decode
	unsigned error = lodepng::load_file(png, filename);
	if (!error) error = lodepng::decode(image, width, height, png);

	//if there's an error, display it
	if (error) std::cout << "decoder error " << error << ": " << lodepng_error_text(error) << std::endl;
	else std::cout << "Image " << filename << " loaded successfully ! " << std::endl;

	//the pixels are now in the vector "image", 4 bytes per pixel, ordered RGBARGBA..., use it as texture, draw it, ...
	png.clear();
	png.shrink_to_fit();
	return;
}

//Encode from raw pixels to an in-memory PNG file first, then write it to disk
//The image argument has width * height RGBA pixels or width * height * 4 bytes
void writeImage(vector<std::vector<unsigned char> >& image, const char* filename, bool color)
{
	std::vector<unsigned char> png, tmp;
	int height = image.size();
	int width = image[0].size();

	// Linearization

	if (!color) // If the input image is a gray image
	{
		for (int j = 0; j < height; j++)
		{
			for (int i = 0; i < width; i++)
			{
				tmp.push_back(image[j][i]);			 // R
				tmp.push_back(image[j][i]);			 // G
				tmp.push_back(image[j][i]);			 // B
				tmp.push_back((unsigned char)255); // A
			}
		}
	}
	else // In case the input image is a colored image with RGBA
	{
		for (int j = 0; j < height; j++)
		{
			for (int i = 0; i < width; i++)
			{
				tmp.push_back(image[j][i]);
			}
		}
	}

	unsigned error = lodepng::encode(png, tmp, width, height);
	if (!error) lodepng::save_file(png, filename);

	//if there's an error, display it
	if (error) std::cout << "encoder error " << error << ": " << lodepng_error_text(error) << std::endl;
	else std::cout << "Image was successfully written to disk at: " << filename << std::endl;

	png.clear();
	tmp.clear();
	png.shrink_to_fit();
	tmp.shrink_to_fit();
	return;
}


// Transform image to gray and downsample it
// It can be applied either on the original or the downsampled image
void rgba2gray(vector<unsigned char>& image, vector<vector<unsigned char>>& gray2D, unsigned& width, unsigned& height, int M)
{
	const unsigned char RGBA = 4;
	M = (int)sqrt(M); // Downsampling factor

	double g; // represents a single gray pixel value 

	// Since the image is reshaped (stored in image2D vector), 
	// let's start RGB2GREY transformation with downsampling embedded too

	std::vector<unsigned char> tmp; // temporary buffer
	unsigned int hgap = width * RGBA * M;

	// Reshaping from 1D to 2D + calculating gray
	for (unsigned int i = 0; i < height / 4; i++)
	{
		for (unsigned int j = i * hgap; j < (i*hgap) + (width*RGBA); j += (M*RGBA))
		{
			g = ((float)image[j] * 0.2126) + ((float)image[j + 1] * 0.7152) + ((float)image[j + 2] * 0.0722);	// Gray calculation
			tmp.push_back((unsigned char)g);
		}
		gray2D.push_back(tmp);
		tmp.clear();
	}

	height /= M;
	width /= M;
	// Freeing up memory
	image.clear();
	image.shrink_to_fit();
	return;
}
