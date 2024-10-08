/**
  @file
  @brief Terrain attributes that can only be calculated with Tarboton's D-infinity flow metric
  @author Richard Barnes (rbarnes@umn.edu), 2015

  This file implements the D-infinite flow routing method originally described by
  Tarboton (1997). It incorporates minor alterations and additional safe-guards
  described in Barnes (2013, TODO).
*/

#pragma once

#include <richdem/common/Array2D.hpp>
#include <richdem/common/constants.hpp>
#include <richdem/common/ProgressBar.hpp>

#include <cassert>
#include <cmath>
#include <queue>
#include <stdexcept>

namespace richdem {

/**
  @brief  Calculate the surface of a digital elevation model
  @author Jenness (2004), Richard Barnes (rbarnes@umn.edu)

    Calculates the surface area of a digital elevation model by connecting the
    central points of cells with triangles and then calculating the area of the
    portion of each triangle which falls within the focal cell. The method is
    described in detail in Jenness (2004)
    <doi:10.2193/0091-7648(2004)032[0829:CLSAFD]2.0.CO;2>

  @param[in]  &elevations   A grid of elevations
  @param[in]   zscale       DEM is scaled by this factor prior to calculation

  @return The surface area of the digital elevation model
*/
template <class T>
double dem_surface_area(
  const Array2D<T> &elevations,
  const double zscale
){
  ProgressBar progress;

  std::cerr<<"\nA DEM Surface Elevation"<<std::endl;
  std::cerr<<"C Jenness, J.S., 2004. Calculating landscape surface area from digital elevation models. Wildlife Society Bulletin 32, 829--839. doi:10.2193/0091-7648(2004)032[0829:CLSAFD]2.0.CO;2"<<std::endl;

  const auto euc_dist = [](const double a, const double b){ return std::sqrt(std::pow((double)a,2.0)+std::pow((double)b,2.0)); };

  //If calculate cell area is lower than actual area, but greater than "AA minus
  //fudge_factor", we clamp the calculated area to the actual area without
  //raising an error. If the calculated area is lower than the fudge_factor, we
  //raise an alarm.
  const double fudge_factor = 1e-4;

  //Using double as an accumulator here is important! Testing this algorithm
  //using the Boost Numeric Interval library should data such as follows:
  //Single-precision sum            = 2.14851e+09
  //Double-precision sum            = 1.61629e+10
  //Single-precision interval Width = 1.09655e+14
  //Single-precision interval Lower = 1.07436e+09
  //Single-precision interval Upper = 1.09656e+14
  //Double-precision interval Width = 181.906
  //Double-precision interval Lower = 1.61629e+10
  //Double-precision interval Upper = 1.61629e+10

  //The upshot is that there is significant uncertainty associated with the
  //floating-point accumulator while the double accumulator has negligible
  //uncertainty.
  double area = 0;

  const double xdist            = elevations.getCellLengthX();
  const double ydist            = elevations.getCellLengthY();
  const double planar_diag_dist = euc_dist(xdist,ydist);

  //TODO: This algorithm would benefit from GPU integration
  progress.start(elevations.size());
  #pragma omp parallel for reduction(+:area)
  for(int y=0;y<elevations.height();y++){
    progress.update( y*elevations.width() );
    for(int x=0;x<elevations.width();x++){
      if(elevations.isNoData(x,y))
        continue;

      //We sum into `cell_area` rather than `area` so that our values are larger
      //when we add to `area`. This helps prevent small numbers from being
      //"swallowed" by large numbers, and other floating-point stuff.
      double cell_area = 0;

      //Loop through neighbours
      for(int n=1;n<=8;n++){
        //This is the next neighbour, which forms part of the triangle
        int nn = n+1;
        if(nn==9) //Wrap around
          nn = 1;

        //In each triangle one neighbour is in the diagonal direction and one is
        //in a straight direction.
        int dn  = n;   //Diagonal Neighbour
        int ndn = nn;  //Not-diagonal neighbour

        //As we walk around half of the time we'll misidentify the diagonal
        //neighbour, but we can fix that here by swapping the labels we just
        //gave the neighbours
        if(!n8_diag[dn])
          std::swap(dn,ndn);

        const double my_elev = zscale*elevations(x,y);

        //Deal with the possibility that the neighbouring cells do not exist. In
        //this case, we pretend that they do exist and are at the same height as
        //the focal cell.
        double dn_elev;
        if(elevations.inGrid(x+d8x[dn],y+d8y[dn]) && !elevations.isNoData(x+d8x[dn],y+d8y[dn]))
          dn_elev = zscale*elevations(x+d8x[dn],y+d8y[dn]);
        else
          dn_elev = my_elev;

        //Do the same for the other neighbour
        double ndn_elev;
        if(elevations.inGrid(x+d8x[ndn],y+d8y[ndn]) && !elevations.isNoData(x+d8x[ndn],y+d8y[ndn]))
          ndn_elev = zscale*elevations(x+d8x[ndn],y+d8y[ndn]);
        else
          ndn_elev = my_elev;

        const double planar_dist_dn   = planar_diag_dist;             //Distance focal cell to diagonal neighbour
        const double planar_dist_ndn  = (d8y[ndn] == 0)?xdist:ydist;  //Distance focal cell to non-diagonal neighbour
        const double planar_dist_bn   = (d8y[ndn] == 0)?ydist:xdist;  //Distance between the neighbour cells

        const double elev_diff_dn     = dn_elev -my_elev; //Elevation drop between focal and diagonal
        const double elev_diff_ndn    = ndn_elev-my_elev; //Elevation drop between focal and non-diagonal
        const double elev_diff_bn     = ndn_elev-dn_elev; //Elevation drop between neighbours

        //Divide these distances by two to form a similar triangle constrained
        //by the boundary of the focal cell
        const double surf_dist_dn     = euc_dist(planar_dist_dn,elev_diff_dn)/2;   //3-space distance between center of focal and diagonal neighbour
        const double surf_dist_ndn    = euc_dist(planar_dist_ndn,elev_diff_ndn)/2; //3-space distance between center of focal and non-diagonal neighbour
        const double surf_dist_bn     = euc_dist(planar_dist_bn,elev_diff_bn)/2;   //3-space distance between neighbours

        //Used to get area of triangle
        const double s = (surf_dist_dn+surf_dist_ndn+surf_dist_bn)/2;

        //Accumulate area of triangle to
        const double tri_area = std::sqrt(s*(s-surf_dist_dn)*(s-surf_dist_ndn)*(s-surf_dist_bn));

        cell_area += tri_area;
      }

      if(cell_area<elevations.getCellArea()){
        if(cell_area+fudge_factor>=elevations.getCellArea())
          cell_area = elevations.getCellArea();
        else
          throw std::runtime_error("A cell had a topographic surface area less than its planar surface area!");
      }

      area += cell_area;
    }
  }
  std::cerr<<"p Succeeded in = "<<progress.stop()<<" s"<<std::endl;

  const double dem_planar_area = elevations.numDataCells()*elevations.getCellArea();
  if(area<dem_planar_area){
    std::cerr<<"W Topographic surface area ("+std::to_string(area)+") < planar surface area ("+std::to_string(dem_planar_area) +")! Choosing planar area.";
    return dem_planar_area;
  }

  return area;
}



enum class PerimType {
  CELL_COUNT,    ///< Counts # of cells bordering DEM edges or NoData cells
  SQUARE_EDGE,   ///< Adds all cell edges bordering DEM edges or NoData cells
};

/**
  @brief  Calculates the perimeter of a digital elevation model
  @author Richard Barnes (rbarnes@umn.edu)

  @param[in]  &arr
  @param[in]  perim_type  A `PerimType` value indicating how to calculate the
                          perimeter.

  @return The perimeter of the digital elevation model
*/
template <class T>
double Perimeter(
  const Array2D<T> &arr,
  const PerimType perim_type
){
  ProgressBar progress;

  std::cerr<<"\nA DEM Perimeter"<<std::endl;
  std::cerr<<"C TODO"<<std::endl;

  unsigned int vertical_edges   = 0;
  unsigned int horizontal_edges = 0;
  unsigned int cell_edges       = 0;

  //TODO: This algorithm would benefit from GPU integration
  progress.start(arr.size());
  #pragma omp parallel for reduction(+:cell_edges) reduction(+:vertical_edges) reduction(+:horizontal_edges)
  for(int y=0;y<arr.height();y++){
    progress.update( y*arr.width() );
    for(int x=0;x<arr.width();x++){
      if(arr.isNoData(x,y))
        continue;

      if(perim_type==PerimType::CELL_COUNT){
        for(int n=1;n<=8;n++){
          if(!arr.inGrid(x+d8x[n],y+d8y[n])){
            cell_edges++;
            break;
          }
        }
      } else if(perim_type==PerimType::SQUARE_EDGE){
        for(int n=1;n<=8;n++){
          if(!arr.inGrid(x+d8x[n],y+d8y[n]) || arr.isNoData(x+d8x[n],y+d8y[n])){
            if(d8x[n]==0) //Pointing at a cell above or below, so horizontal edge
              horizontal_edges++;
            else if(d8y[n]==0) //Point at cell left or right, so vertical edge
              vertical_edges++;
          }
        }
      } else {
        throw std::runtime_error("Unrecognised PerimType!");
      }
    }
  }
  std::cerr<<"p Succeeded in = "<<progress.stop()<<" s"<<std::endl;

  if(perim_type==PerimType::CELL_COUNT)
    return cell_edges;
  else if(perim_type==PerimType::SQUARE_EDGE)
    return horizontal_edges*arr.getCellLengthX()+vertical_edges*arr.getCellLengthY();
  else
    throw std::runtime_error("Unrecognised PerimType!");
}



/**
  @brief  Applies a bucket-fill paint operation to one raster based on another
  @author Richard Barnes (rbarnes@umn.edu)

  @param[in]      &check_raster  Raster whose values are checked for the BucketFill
  @param[in,out]  &set_raster    Raster whose values are set by the BucketFill.
                                 If `set_raster` already has `set_value`, then
                                 the FloodFill won't progress over it. This
                                 avoids needing a separate visisted raster.
  @param[in]      check_value  Value in `check_raster` which indicates a value
                               in `set_raster` should be set
  @param[in]      set_value    Value that `set_raster` is set to
  @param[in]      &seed        Vector of seed cells to seed the BucketFill
*/
template<Topology topo, class T, class U>
void BucketFill(
  const Array2D<T>    &check_raster,
  Array2D<U>          &set_raster,
  const T             &check_value,
  const U             &set_value,
  std::vector<size_t> &seeds
){
  if(check_raster.width()!=set_raster.width() || check_raster.height()!=set_raster.height()){
    throw std::runtime_error("Rasters must have the same dimension for BucketFill!");
  }

  static_assert(topo==Topology::D8 || topo==Topology::D4);
  constexpr auto dx = get_dx_for_topology<topo>();
  constexpr auto dy = get_dy_for_topology<topo>();
  constexpr auto nmax = get_nmax_for_topology<topo>();

  while(!seeds.empty()){
    const auto c = seeds.back();
    seeds.pop_back();

    if(check_raster(c)!=check_value || set_raster(c)==set_value){
      continue;
    }

    set_raster(c) = set_value;

    const auto [cx, cy] = check_raster.iToxy(c);

    for(int n=1;n<=nmax;n++){
      if(!check_raster.inGrid(cx+dx[n],cy+dy[n]))
        continue;

      const auto ni = check_raster.xyToI(cx+dx[n], cy+dy[n]);

      if(check_raster(ni)==check_value && set_raster(ni)!=set_value)
        seeds.emplace_back(ni);
    }
  }
}



/**
  @brief  Applies a bucket-fill paint operation to one raster based on another
          starting from the edges.
  @author Richard Barnes (rbarnes@umn.edu)

  @param[in]      &check_raster  Raster whose values are checked for the BucketFill
  @param[in,out]  &set_raster    Raster whose values are set by the BucketFill.
                                 If `set_raster` already has `set_value`, then
                                 the FloodFill won't progress over it. This
                                 avoids needing a separate visisted raster.
  @param[in]      check_value  Value in `check_raster` which indicates a value
                               in `set_raster` should be set
  @param[in]      set_value    Value that `set_raster` is set to
*/
template<Topology topo, class T, class U>
void BucketFillFromEdges(
  const Array2D<T>   &check_raster,
  Array2D<U>         &set_raster,
  const T            &check_value,
  const U            &set_value
){
  std::vector<size_t> seeds;
  seeds.reserve(2*check_raster.width()+2*check_raster.height());
  for(int y=0;y<check_raster.height();y++){
    seeds.emplace_back(check_raster.xyToI(0,y));
    seeds.emplace_back(check_raster.xyToI(check_raster.width()-1,y));
  }
  for(int x=0;x<check_raster.width();x++){
    seeds.emplace_back(check_raster.xyToI(x,0));
    seeds.emplace_back(check_raster.xyToI(x,check_raster.height()-1));
  }
  BucketFill<topo>(check_raster, set_raster, check_value, set_value, seeds);
}

}
