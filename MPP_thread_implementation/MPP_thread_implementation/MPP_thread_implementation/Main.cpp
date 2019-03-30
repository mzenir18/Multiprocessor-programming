#include <iostream>
#include <vector>
#include <fstream>
#include <numeric>
#include <chrono>
#include <fstream> // To be removed later
#include <omp.h> // OpenPM for multi-threading
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
int WIN_SIZE = 9; // Window size
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
inline void depthEstim(img_manipulator& img1, img_manipulator& img2);
void crossCheck(img_manipulator& img1, img_manipulator& img2);
void occlusionFilling(vector<vector <unsigned char> >& image);
/*
	########################################## Main function ##########################################
*/

int main()
{
	int list_nbr_threads[] = { 1, 2 , 3, 4, 6, 8, 12, 16, 32, 512, 1024, 1536, 2048, 350000 };

	int list_win_size[] = { 9, 13, 17, 21};

	// Output results file 
	ofstream resultsFile;
	resultsFile.open("TH_results/Inner_Loop_1.txt"); // Please change the name to avoid overwriting current files
	
	for (int th : list_nbr_threads)
	{
		// Set thread counts
		N = th;
		for (int w: list_win_size)
		{
			// Set window size
			WIN_SIZE = w;

			// Calculating execution time
			auto start = chrono::steady_clock::now();

			// // Loading images to memory
			// Creating an image manipulator structure object
			img_manipulator img1;
			img_manipulator img2;

			// Let's load the images into memory
			loadImages(img1.image, img1.width, img1.height, FILE_1); // Loading first image
			loadImages(img2.image, img2.width, img2.height, FILE_2); // Loading second image

			// Print out inital loaded image characteristics
			///cout << "Height of img1: " << img1.height << " Width of img1: " << img1.width << endl;
			///cout << "Image: " << img1.image.size() << endl;

			// Gray2D
			///cout << "Applying gray on image1... " << endl;
			rgba2gray(img1.image, img1.gray2D, img1.width, img1.height, M);
			///cout << "Applying gray on image2... " << endl;
			rgba2gray(img2.image, img2.gray2D, img2.width, img2.height, M);

			// Print out inital loaded image characteristics
			///cout << "Height of gray2D 1 " << img1.gray2D.size() << " Width of img1: " << img1.gray2D[0].size() << endl;
			///cout << "Height of gray2D 2  " << img2.gray2D.size() << " Width of img2: " << img2.gray2D[0].size() << endl;

			// Applying depth estimation using ZNCC
	
			///cout << "Executing depth estimation function ..." << endl;
			depthEstim(img1, img2);
	
			// Writting disparity maps to disk
			writeImage(img1.disp_img, FILE_L2R, false); // (2D image, filename, IsItColored)
			// Writting disparity maps to disk
			writeImage(img2.disp_img, FILE_R2L, false); // (2D image, filename, IsItColored)

			// Cross check
			///cout << "Cross check running ..." << endl;
			crossCheck(img1, img2);
			// Writting cross check output to disk
			writeImage(img1.disp_img, FILE_CROSS_CHECKED, false);
	
			// Occlusion filling
			///cout << "Occlusion filling running ..." << endl;
			occlusionFilling(img1.disp_img);
			// Writting occlusion filled final image to disk
	
			writeImage(final_img, FILE_OCCLUSION_FILLED, false); // (2D image, filename, IsItColored)
	
			// End of program time measurement
			auto end = chrono::steady_clock::now();

			// Store the time difference between start and end
			auto diff = end - start;

			// Printing out the execution duration
			///cout << chrono::duration <double, milli>(diff).count() / 1000 << " s" << endl;
			cout << "Thread count: " << N << " ; WinSize: " << WIN_SIZE << " ; Execution time: " << chrono::duration <double, milli>(diff).count() / 1000 << " s" << endl;
			// Writing data to file
			resultsFile << "Thread count: " << N << " ; WinSize: " << WIN_SIZE << " ; Execution time: " << chrono::duration <double, milli>(diff).count() / 1000 << " s\n";
		}
	}
	resultsFile.close();
	system("pause");
	return 0;
}

/*
	---------------------------------------- Function implementations -------------------------------------------------
*/

// Source: https://raw.githubusercontent.com/lvandeve/lodepng/master/examples/example_decode.cpp
//Load PNG file from disk to memory first, then decode to raw pixels in memory.
void loadImages(vector<unsigned char>& image, unsigned& width, unsigned& height, const char* filename) {
	std::vector<unsigned char> png; // Image will be represented in 8-bit unsinged representation

	//load and decode
	unsigned error = lodepng::load_file(png, filename);
	if (!error) error = lodepng::decode(image, width, height, png);

	//if there's an error, display it
	///if (error) std::cout << "decoder error " << error << ": " << lodepng_error_text(error) << std::endl;
	///else std::cout << "Image " << filename << " loaded successfully ! " << std::endl;

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
	///if (error) std::cout << "encoder error " << error << ": " << lodepng_error_text(error) << std::endl;
	///else std::cout << "Image was successfully written to disk at: " << filename << std::endl;

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


// Depth estimation impelmentation
// Parallel execution
inline void depthEstim(img_manipulator& img1, img_manipulator& img2)
{
	unsigned int best_disp1{ 0 }, best_disp2{ 0 };
	int xs1{ 0 }, xe1{ 0 }, xs2{ 0 }, xe2{ 0 };
	double window_mean_1{ 0 }, window_mean_2{ 0 }, zncc_value1{ 0 }, current_max_sum1{ 0 }, zncc_value2{ 0 }, current_max_sum2{ 0 };
	int edge = (int)(WIN_SIZE / 2);
	const unsigned MAPPING_FACTOR = (256 / MAX_DISP);
	int display{ 50 };

	vector<vector<unsigned char>> disp_map1(img1.height - edge*2, vector<unsigned char>(img1.width - edge*2, (unsigned char)0)); // Initialize disp map to zero for later indexing
	vector<vector<unsigned char>> disp_map2(img2.height - edge*2, vector<unsigned char>(img2.width - edge*2, (unsigned char)0));

	omp_set_num_threads(N);
	// Outer loop parallelization
#pragma omp parallel for default(none) firstprivate(best_disp1, best_disp2, xs1, xe1, xs2, xe2, window_mean_1, window_mean_2, zncc_value1, zncc_value2, current_max_sum2, edge, MAPPING_FACTOR, display, WIN_SIZE) shared(img1, img2, disp_map1, disp_map2)
	for (int j = edge; j < img1.height - edge; j++)
	{
		/*
		if ((j % display) == 0)
		{
			cout << "j: " << j << endl;
			printf("Thread: %d\n", omp_get_thread_num());
		}
		*/
// Uncomment line below to parallelize the inner for loop
//#pragma omp parallel for default(none) firstprivate(best_disp1, best_disp2, xs1, xe1, xs2, xe2, window_mean_1, window_mean_2, zncc_value1, zncc_value2, current_max_sum2, edge, MAPPING_FACTOR, display, WIN_SIZE) shared(j, img1, img2, disp_map1, disp_map2)
// Uncomment line below to implement both parallel for loops
//#pragma omp parallel for default(none) shared(best_disp1, best_disp2, xs1, xe1, xs2, xe2, window_mean_1, window_mean_2, zncc_value1, zncc_value2, current_max_sum2, edge, MAPPING_FACTOR, display, WIN_SIZE, j, img1, img2, disp_map1, disp_map2)
		for (int i = edge; i < img1.width - edge; i++)
		{

			//######## Mean calculation starts here
			window_mean_1 = 0;
			window_mean_2 = 0;
			// Window mean calculation
			for (int k = -edge; k < edge + 1; k++)
			{
				for (int m = -edge; m < edge + 1; m++)
				{
					window_mean_1 += (double)img1.gray2D[j + k][i + m];
					window_mean_2 += (double)img2.gray2D[j + k][i + m];
				}
			}
			window_mean_1 /= (WIN_SIZE*WIN_SIZE); // Mean of the window from image 1
			window_mean_2 /= (WIN_SIZE*WIN_SIZE); // Mean of the window from image 2

			//######## Mean calculation ends here !

			current_max_sum1 = -1;
			current_max_sum2 = -1;

			//######## ZNCC calculation starts here !
			for (unsigned d = 0; d < MAX_DISP; d++)
			{
				double nominator1{ 0 }, denominator1{ 0 }, left1{ 0 }, right1{ 0 }, left_den1{ 0 }, right_den1{ 0 };
				double nominator2{ 0 }, denominator2{ 0 }, right2{ 0 }, left2{ 0 }, left_den2{ 0 }, right_den2{ 0 };

				// Checking if the window is inside the image or not
				xs1 = i - edge - d;
				xe1 = i + edge + 1 - d;

				xs2 = i - edge + d;
				xe2 = i + edge + 1 + d;

				if ((xs1 < 0 || xe1 > img1.width - edge) and (xs2 < 0 || xe2 > img2.width - edge)) // Both out of bounds
				{	// For this value of d, window is out of the image
					zncc_value1 = -1;
					zncc_value2 = -1;
				}
				else if (xs1 < 0 || xe1 > img1.width - edge) // Only img1:d is out of bound
				{

					for (int k = -edge; k < edge + 1; k++)
					{
						for (int m = -edge; m < edge + 1; m++)
						{
							left1 = ((double)img2.gray2D[j + k][i + m] - window_mean_2);			// for img2
							right2 = ((double)img1.gray2D[j + k][i + m + d] - window_mean_1);		// for img2
							nominator2 += left1 * right2;
							left_den1 += pow(left1, 2);
							right_den2 += pow(right2, 2);
						}
					}
					denominator2 = sqrt(left_den1) * sqrt(right_den2);

					// zncc value
					zncc_value1 = -1; // Because img1:d window is out of bound
					zncc_value2 = nominator2 / denominator2;

				}
				else if (xs2 < 0 || xe2 > img2.width - edge) // Only img2:d is out of bound
				{
					for (int k = -edge; k < edge + 1; k++)
					{
						for (int m = -edge; m < edge + 1; m++)
						{
							left1 = ((double)img1.gray2D[j + k][i + m] - window_mean_1); // img1
							right1 = ((double)img2.gray2D[j + k][i + m - d] - window_mean_2);		// img1
							nominator1 += left1 * right1;
							left_den1 += pow(left1, 2);
							right_den1 += pow(right1, 2);
						}
					}
					denominator1 = sqrt(left_den1) * sqrt(right_den1);

					// zncc value
					zncc_value1 = nominator1 / denominator1;
					zncc_value2 = -1; // Because img2:d window is out of bound

				}
				else // Both in bound
				{
					for (int k = -edge; k < edge + 1; k++)
					{
						for (int m = -edge; m < edge + 1; m++)
						{
							left1 = ((double)img1.gray2D[j + k][i + m] - window_mean_1);			// img1
							left2 = ((double)img2.gray2D[j + k][i + m] - window_mean_2);			// img2
							right1 = ((double)img2.gray2D[j + k][i + m - d] - window_mean_2);		// img1
							right2 = ((double)img1.gray2D[j + k][i + m + d] - window_mean_1);		// img2
							nominator1 += left1 * right1;
							nominator2 += left2 * right2;
							left_den1 += pow(left1, 2);
							left_den2 += pow(left2, 2);
							right_den1 += pow(right1, 2);
							right_den2 += pow(right2, 2);

						}
					}
					denominator1 = sqrt(left_den1) * sqrt(right_den1);
					denominator2 = sqrt(left_den2) * sqrt(right_den2);

					// zncc value
					zncc_value1 = nominator1 / denominator1;
					zncc_value2 = nominator2 / denominator2;

				}


				//######## ZNCC calculation ends here !

				if (zncc_value1 > current_max_sum1)
				{
					// Updating current maximum sum
					current_max_sum1 = zncc_value1;
					// Updating best disparity value
					best_disp1 = d * MAPPING_FACTOR; // [0,256)

				} // End if

				if (zncc_value2 > current_max_sum2)
				{
					// Updating current maximum sum
					current_max_sum2 = zncc_value2;
					// Updating best disparity value
					best_disp2 = d * MAPPING_FACTOR; // [0,256)

				} // End if


			} // End disp for
			// Adding the elements
			disp_map1[j - edge][i - edge] = best_disp1;
			disp_map2[j - edge][i - edge] = best_disp2;

		} // End width for
	} // End height for

	img1.disp_img = disp_map1;
	img2.disp_img = disp_map2;

	disp_map1.clear();
	disp_map2.clear();
	img1.gray2D.clear();
	img2.gray2D.clear();
	img1.image.clear();
	img2.image.clear();
	return;
}


// Cross check implementation
void crossCheck(img_manipulator& img1, img_manipulator& img2)
{
	// Results will be stored in gray2D
	int abs_diff{ 0 };
	const int edge = (int)((WIN_SIZE / 2)) * 2;

	for (unsigned j = 0; j < img1.height - edge; j++)
	{
		for (unsigned i = 0; i < img1.width - edge; i++)
		{
			abs_diff = abs((int)img1.disp_img[j][i] - (int)img2.disp_img[j][i]);
			if (abs_diff > THRESHOLD) // If difference is too high
			{
				if ((int)img1.disp_img[j][i] > (int)img2.disp_img[j][i])
					img1.disp_img[j][i] = img1.disp_img[j][i];
				else
					img1.disp_img[j][i] = img2.disp_img[j][i];
			}
			else // Small difference, we will change the pixel by the mean
			{
				img1.disp_img[j][i] = (img1.disp_img[j][i] + img2.disp_img[j][i]) / 2;
			}
		}
	}
	img2.disp_img.clear();
	// Updating the height and the width
	img1.height = img1.disp_img.size();
	img1.width = img1.disp_img[0].size();
	return;
}


// Occlusion filling
void occlusionFilling(vector<vector <unsigned char> >& image)
{
	int height = image.size();
	int width = image[0].size();
	int direction{ -1 }, counter{ 0 }, tmp;

	for (int j = 0; j < height; j++)
	{
		for (int i = 0; i < width; i++)
		{
			if ((int)image[j][i] < 20) // Pixel with zero => let's fill it
			{
				tmp = i; // Saving the value of i
				while (true)
				{
					tmp += direction; // Going left (assigning the closest non zero value left neighbor) (direction is negative) 

					if (tmp < 0) // Check that we are going in the correct direction
					{
						// No neighbor on the left, we should go write
						tmp = i;
						tmp++;
						direction = 1;
					}
					else if (tmp > width)
					{
						cout << " i : " << i << " j : " << j << endl;
						system("pause");
						break; // It shoud never reach this statement
					}

					if ((int)image[j][tmp] > 20)
					{
						image[j][i] = image[j][tmp];  // Assigning new value to the zero pixel
						direction = -1; // Reset the direction
						break;
					}
				}
			}
		}

	}
	final_img = image;
	image.clear();
	return;
}