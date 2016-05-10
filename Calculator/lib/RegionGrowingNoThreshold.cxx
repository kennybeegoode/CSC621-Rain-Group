
/*=========================================================================
 *
 *  Copyright Insight Software Consortium
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/
#include "RegionGrowingNoThreshold.h"

#include <vector>
#include <unordered_map>

#include "itkConnectedThresholdImageFilter.h"
#include "itkImage.h"
#include "itkCastImageFilter.h"
#include "itkRGBPixel.h"
#include "itkLaplacianRecursiveGaussianImageFilter.h"
#include "itkSubtractImageFilter.h"
#include "itkCurvatureAnisotropicDiffusionImageFilter.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#define MAX(a, b) ((a) >= (b) ? (a) : (b))
#define MIN(a, b) ((a) <= (b) ? (a) : (b))

// Compute statistics variables.
void RegionGrowingNoThreshold::ComputeStats(InternalImageType::Pointer& image,
                  std::unordered_map<long long, InternalImageType::IndexType>& region_points,
                  double* mgv, double* upper_dev, double* lower_dev) {
  long long sum = 0;
  // Computing mean for all points
  for (const auto& kv : region_points) {
    sum += image->GetPixel(kv.second);
  }
  *mgv = sum / region_points.size();

  // Computing mean for upper and lower halves.
  long long upper_sum = 0, lower_sum = 0;
  long long upper_num = 0, lower_num = 0;
  for (const auto& kv : region_points) {
    if (image->GetPixel(kv.second) >= *mgv) {
      upper_sum += image->GetPixel(kv.second);
      upper_num++;
    }
    if (image->GetPixel(kv.second) <= *mgv) {
      lower_sum += image->GetPixel(kv.second);
      lower_num++;
    }
  }
  double upper_mean = upper_sum / upper_num;
  double lower_mean = lower_sum / lower_num;

  // Computing standard deviation for upper and lower halves
  *upper_dev = 0;
  *lower_dev = 0;
  for (const auto& kv : region_points) {
    if (image->GetPixel(kv.second) >= *mgv) {
      *upper_dev += (image->GetPixel(kv.second) - upper_mean) * (image->GetPixel(kv.second) - upper_mean);
    }
    if (image->GetPixel(kv.second) <= *mgv) {
      *lower_dev += (image->GetPixel(kv.second) - lower_mean) * (image->GetPixel(kv.second) - lower_mean);
    }
  }

  *upper_dev /= upper_num;
  *lower_dev /= lower_num;

  *upper_dev = sqrt(*upper_dev);
  *lower_dev = sqrt(*lower_dev);
}

// Compute thresholds for the first run.
void RegionGrowingNoThreshold::ComputeThreshold(
                      InternalImageType::Pointer& image,
                      std::unordered_map<long long, InternalImageType::IndexType>& region_points,
                      double* upper_threshold, double* lower_threshold) {
  double upper_dev, lower_dev, mgv;
  ComputeStats(image, region_points, &mgv, &upper_dev, &lower_dev);

  *upper_threshold = mgv + (upper_dev * 1.5 + 20 / sqrt(region_points.size()));
  *lower_threshold = mgv - (lower_dev * 1.5 + 20 / sqrt(region_points.size()));

  //std::cout << "mgv " << mgv << " upper_dev " << upper_dev << " lower_dev " << lower_dev << std::endl;
  //std::cout << "thresholds " << *upper_threshold << " " << *lower_threshold << std::endl;
}

// Compute thresholds for the final run.
void RegionGrowingNoThreshold::ComputeFinalThreshold(
                           InternalImageType::Pointer& image,
                           std::unordered_map<long long, InternalImageType::IndexType>& region_points,
                           double* upper_threshold, double* lower_threshold, double* upper_dev,
                           double* lower_dev) {
  double mgv;
  ComputeStats(image, region_points, &mgv, upper_dev, lower_dev);

  if (*lower_dev == 0) *lower_dev = 1;
  if (*upper_dev == 0) *upper_dev = 1;

  *upper_threshold = mgv + (*upper_dev * 2.3 * 2.58 + 20 / sqrt(region_points.size()));
  *lower_threshold = mgv - (*lower_dev * 2.3 * 2.58 + 20 / sqrt(region_points.size()));

  std::cout << "mgv " << mgv << " upper_dev " << *upper_dev << " lower_dev " << *lower_dev << std::endl;
  std::cout << "final thresholds " << *upper_threshold << " " << *lower_threshold << std::endl;
}

// Compute unique value for every point.
long long RegionGrowingNoThreshold::ComputePointHash(InternalImageType::IndexType& point) {
  return point[0] * point[1] * point[2] + point[0]*100 + point[1]*10 + point[2];
}

// Check if point is within bounds.
bool RegionGrowingNoThreshold::CheckBounds(
                 InternalImageType::IndexType& point,
                 InternalImageType::SizeType size) {
  if (!(point[0] < size[0] && point[0] >= 0)) return false;
  if (!(point[1] < size[1] && point[1] >= 0)) return false;
  if (!(point[2] < size[2] && point[2] >= 0)) return false;
  return true;
}

bool RegionGrowingNoThreshold::GrowRegions(
                 InternalImageType::Pointer& image,
                 InternalImageType::IndexType& seed,
                 double* upper_threshold,
                 double* lower_threshold) {
  std::unordered_map<long long, InternalImageType::IndexType> visited_points;
  std::queue<InternalImageType::IndexType> q;
  InternalImageType::RegionType region = image->GetLargestPossibleRegion();
  InternalImageType::SizeType size = region.GetSize();

  // Adding the first point to the region.
  visited_points[ComputePointHash(seed)] = seed;
  // Adding all 26 neighbors of the seed point
  for (int i = -1; i <= 1; i++) {
    for (int j = -1; j <= 1; j++) {
      for (int k = -1; k <= 1; k++) {
        // We should not add cental point
        if (!(i == 0 && j == 0 && k ==0)) {
          InternalImageType::IndexType seed_copy = seed;
          seed_copy[0] += i;
          seed_copy[1] += j;
          seed_copy[2] += k;
          if (!CheckBounds(seed_copy, size)) continue;
          q.push(seed_copy);
          //std::cout << "new point " << image->GetPixel(seed_copy) << std::endl;
          visited_points[ComputePointHash(seed_copy)] = seed_copy;
        }
      }
    }
  }
  ComputeThreshold(image, visited_points, upper_threshold, lower_threshold);
  long long region_size = visited_points.size();

  // Running a loop - removing front element from the queue and adding its neighbors
  // to the back of the queue.
  int iter = 0;
  while (q.size()) {
    iter++;
    //std::cout << "iter " << iter << " " << visited_points.size() << " queue " << q.size() << std::endl;
    if (iter > 10000000) {
      break;
    }
    InternalImageType::IndexType elem = q.front();
    q.pop();
    for (int i = -1; i <= 1; i++) {
      for (int j = -1; j <= 1; j++) {
        for (int k = -1; k <= 1; k++) {
          // We should not add cental point
          if (!(i == 0 && j == 0 && k ==0)) {
            InternalImageType::IndexType seed_copy = elem;
            seed_copy[0] += i;
            seed_copy[1] += j;
            seed_copy[2] += k;
            if (!CheckBounds(seed_copy, size)) continue;
            // If pixel value is withing the thresholds, we proceed with it.
            //std::cout << "seed_copy " << seed_copy << " " << image->GetPixel(seed_copy) << std::endl;

            if (image->GetPixel(seed_copy) >= *lower_threshold && image->GetPixel(seed_copy) <= *upper_threshold) {
              //std::cout << "fits in thresh" << std::endl;
              // If we already visited this point, don't add it.
              if (visited_points.find(ComputePointHash(seed_copy)) !=
                  visited_points.end()) {
                continue;
              }
              //std::cout << "new point " << image->GetPixel(seed_copy) << std::endl;

              q.push(seed_copy);
              visited_points[ComputePointHash(seed_copy)] = seed_copy;

              // Recomputing thresholds if number of points doubled.
              if (visited_points.size() / region_size >= 2) {
                ComputeThreshold(image, visited_points, upper_threshold, lower_threshold);
                region_size = visited_points.size();
              }
            }
          }
        }
      }
    }
  }
  double upper_dev, lower_dev;
  // Finally, we need to compute the final thresholds for the second run.
  ComputeFinalThreshold(image, visited_points, upper_threshold, lower_threshold,
                        &upper_dev, &lower_dev);
  // If deviation was too big, it returns false.
  if (lower_dev >= 9 || upper_dev >= 9) {
    return false;
  }
  return true;
}

void RegionGrowingNoThreshold::ComputePerSliceAvgAndDev(InternalImageType::Pointer& image,
  double* avg, double* dev) {
  InternalImageType::RegionType region = image->GetLargestPossibleRegion();
  InternalImageType::SizeType size = region.GetSize();
  long long num_slices = 0;
  std::cout << "avg start " << std::endl;
  for (int z = 0; z < size[2]; z++) {
    long long num_points = 0;
    for (int x = 0; x < size[0]; x++) {
      for (int y = 0; y < size[1]; y++) {
        InternalImageType::IndexType pt;
        pt[0] = x;
        pt[1] = y;
        pt[2] = z;
        // Found a point which has pixel.
        if (image->GetPixel(pt) > 0) {
          num_points++;
        }
      }
    }
    if (num_points > 0) {
      *avg += num_points;
      num_slices++;
    }
  }
  if (num_slices > 0) *avg /= num_slices;
  std::cout << "avg end " << *avg << std::endl;

  std::cout << "dev start " << std::endl;

  for (int z = 0; z < size[2]; z++) {
    long long num_points = 0;
    for (int x = 0; x < size[0]; x++) {
      for (int y = 0; y < size[1]; y++) {
        InternalImageType::IndexType pt;
        pt[0] = x;
        pt[1] = y;
        pt[2] = z;
        // Found a point which has pixel.
        if (image->GetPixel(pt) > 0) {
          num_points++;
        }
      }
    }
    if (num_points > 0) {
      *dev += (num_points - *avg) * (num_points - *avg);
      num_slices++;
    }
  }
  if (num_slices > 0) *dev = sqrt(*dev / num_slices);
  std::cout << "dev end " << *dev << std::endl;

}


std::vector<std::vector<int>> RegionGrowingNoThreshold::ComputeCentroids(
  InternalImageType::Pointer& image, double* max_image_x_dist_ratio,
  double* max_image_y_dist_ratio, double* max_image_z_dist_ratio, bool final) {
  std::vector<std::vector<int>> centroids;
  InternalImageType::RegionType region = image->GetLargestPossibleRegion();
  InternalImageType::SizeType size = region.GetSize();
  *max_image_x_dist_ratio = -1;
  *max_image_y_dist_ratio = -1;
  *max_image_z_dist_ratio = -1;
  long long z_max = -1, z_min = 100000;
  double avg = 0, dev;
  if (final) {
    ComputePerSliceAvgAndDev(image, &avg, &dev);
  } else {
    dev = 100000000000;
  }
  // Running over all slices and x an y coordinates.
  for (int z = 0; z < size[2]; z++) {
      long long x_max = -1, x_min = 100000, y_max = -1, y_min = 100000;
      long long x_sum = 0, y_sum = 0;
      long long num = 0;
      for (int x = 0; x < size[0]; x++) {
        for (int y = 0; y < size[1]; y++) {
          InternalImageType::IndexType pt;
          pt[0] = x;
          pt[1] = y;
          pt[2] = z;
          // Found a point which has pixel.
          if (image->GetPixel(pt) > 0) {
            x_max = MAX(x, x_max);
            x_min = MIN(x, x_min);
            y_max = MAX(y, y_max);
            y_min = MIN(y, y_min);
            z_max = MAX(z, z_max);
            z_min = MIN(z, z_min);
            *max_image_x_dist_ratio = MAX(x_max - x_min, *max_image_x_dist_ratio);
            *max_image_y_dist_ratio = MAX(y_max - y_min, *max_image_y_dist_ratio);
            *max_image_z_dist_ratio = MAX(z_max - z_min, *max_image_z_dist_ratio);
            x_sum += x;
            y_sum += y;
            num++;
          }
        }
      }
      // Compute centroid if number of points in the slice was non-zero and it falls
      // within the deviation interval.
      if (num != 0 && num < avg + 2 * dev && num > avg - 2 * dev) {
        std::vector<int> coord;
        coord.push_back(x_sum / num);
        coord.push_back(y_sum / num);
        coord.push_back(z);
        centroids.push_back(coord);
      }
  }

  *max_image_x_dist_ratio = *max_image_x_dist_ratio / size[0];
  *max_image_y_dist_ratio = *max_image_y_dist_ratio / size[1];
  *max_image_z_dist_ratio = *max_image_z_dist_ratio / size[2];
  std::cout << "max_x_dist_ratio " << *max_image_x_dist_ratio
            << " max_y_dist_ratio " << *max_image_y_dist_ratio
            << " max_z_dist_ratio " << *max_image_z_dist_ratio << std::endl;

  return centroids;
}

std::vector<std::vector<int>> RegionGrowingNoThreshold::GetCentroids(char filename[], int seed_x, int seed_y, int seed_z) {
  typedef  itk::ImageFileReader< InternalImageType > ReaderType;
  //typedef itk::CastImageFilter< InternalImageType, OutputImageType >
  //                                                 CastingFilterType;

  ReaderType::Pointer reader = ReaderType::New();
  //CastingFilterType::Pointer caster = CastingFilterType::New();

  reader->SetFileName( filename );

  // Smoothing the image before appplying the region growing.
  typedef itk::CurvatureAnisotropicDiffusionImageFilter< InternalImageType, InternalImageType >
     CurvatureAnisotropicDiffusionImageFilterType;
  CurvatureAnisotropicDiffusionImageFilterType::Pointer smoothing =
     CurvatureAnisotropicDiffusionImageFilterType::New();
  smoothing->SetNumberOfIterations(10 );
  smoothing->SetTimeStep(0.020);
  smoothing->SetConductanceParameter(3);

  typedef itk::LaplacianRecursiveGaussianImageFilter< InternalImageType, InternalImageType >
    LaplacianRecursiveGaussianImageFilterType;
  LaplacianRecursiveGaussianImageFilterType::Pointer laplacian =
    LaplacianRecursiveGaussianImageFilterType::New();

  typedef itk::SubtractImageFilter<InternalImageType> SubtractImageFilterType;
  SubtractImageFilterType::Pointer diff = SubtractImageFilterType::New();

  smoothing->SetInput(reader->GetOutput());
  laplacian->SetInput(smoothing->GetOutput());
  diff->SetInput1(smoothing->GetOutput());
  diff->SetInput2(laplacian->GetOutput());
  reader->Update();
  InternalImageType::Pointer image = diff->GetOutput();

  //typedef  itk::ImageFileWriter<  OutputImageType  > WriterType;
  //WriterType::Pointer writer = WriterType::New();
  //writer->SetFileName("output1.mhd");

  typedef itk::ConnectedThresholdImageFilter< InternalImageType,
       InternalImageType > ConnectedFilterType;
  ConnectedFilterType::Pointer connectedThreshold = ConnectedFilterType::New();
  // The color we are going to paint the region with.
  connectedThreshold->SetReplaceValue(255);

  // Setting up seeds.
  InternalImageType::IndexType  index;
  index[0] = seed_x;
  index[1] = seed_y;
  index[2] = seed_z;

  double max_image_x_dist_ratio, max_image_y_dist_ratio,
      max_image_z_dist_ratio;
  diff->Update();
  long long i_max = 0;
  long long j_max = 0;
  long long k_max = 0;
  double z_max = 0;

  // Going over all neibghorhood if the seed and trying each of them as the
  // seed, at the same time ,we make sure that we did not do oversegmentation
  // and that the length of the spine is not too short.
  for (int i = -2; i <= 2; i++) {
    for (int j = -2; j <= 2; j++) {
      for (int k = -1; k <= 1; k++) {
        InternalImageType::IndexType seed_copy = index;
        seed_copy[0] += i;
        seed_copy[1] += j;
        seed_copy[2] += k;
        connectedThreshold->SetSeed(seed_copy);
        InternalImageType::Pointer image = reader->GetOutput();

        double upper_threshold, lower_threshold;
        // GrowRegions returning false means that deviation of points was too
        // big.
        // The (i == 3 && j == 3 && k == 3) statement makes sure that even
        // when deviation is too big for all neighbors, we still produce some
        // output.
        if (!GrowRegions(image, seed_copy, &upper_threshold, &lower_threshold) &&
            !(i == 2 && j == 2 && k == 1)){
          continue;
        }

        connectedThreshold->SetLower(lower_threshold);
        connectedThreshold->SetUpper(upper_threshold);
        connectedThreshold->SetInput(image);
        connectedThreshold->Update();

        image = connectedThreshold->GetOutput();
        double max_image_x_dist_ratio, max_image_y_dist_ratio,
            max_image_z_dist_ratio;
        ComputeCentroids(image, &max_image_x_dist_ratio, &max_image_y_dist_ratio,
                         &max_image_z_dist_ratio, false);
        // If we don't have oversegmentation, we check if it is the longest
        // segment we ever saw and if it is the case record the shift values
        if ((max_image_x_dist_ratio <= 0.33 && max_image_y_dist_ratio <= 0.33 &&
            max_image_z_dist_ratio >= 0.6) || (i == 2 && j == 2 && k == 1)) {
          std::cout << "break" << std::endl;
          if (max_image_z_dist_ratio > z_max) {
            z_max = max_image_z_dist_ratio;
            i_max = i;
            j_max = j;
            k_max = k;
            i = 2;
            j = 2;
            k = 1;
          }
        }
      }
    }
  }
  // Doing the final segmentation using shift values maximizing the spine
  // height.
  InternalImageType::IndexType seed_copy = index;
  seed_copy[0] += i_max;
  seed_copy[1] += j_max;
  seed_copy[2] += k_max;
  connectedThreshold->SetSeed(seed_copy);
  image = diff->GetOutput();

  double upper_threshold, lower_threshold;
  GrowRegions(image, seed_copy, &upper_threshold, &lower_threshold);

  connectedThreshold->SetLower(lower_threshold);
  connectedThreshold->SetUpper(upper_threshold);
  connectedThreshold->SetInput(image);
  connectedThreshold->Update();

  image = connectedThreshold->GetOutput();

  //caster->SetInput(connectedThreshold->GetOutput());
  //writer->SetInput(caster->GetOutput());

  //writer->Update();

  image = connectedThreshold->GetOutput();
  return ComputeCentroids(image, &max_image_x_dist_ratio, &max_image_y_dist_ratio,
    &max_image_z_dist_ratio, true);
}
