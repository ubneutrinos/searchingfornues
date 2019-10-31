#ifndef LLRPID_H
#define LLRPID_H

#include "LLRPID_proton_muon_lookup.h"

namespace searchingfornues
{

  class LLRPID
  {
  public:

    LLRPID(){}

    void set_dedx_binning(size_t plane, size_t num_bins, std::vector<float> bin_edges)
    {
      dedx_num_bins[plane] = num_bins;
      dedx_bin_edges[plane] = bin_edges;
    }

    void set_par_binning(size_t plane, std::vector<size_t> num_bins, std::vector<std::vector<float>> bin_edges)
    {
      parameters_num_bins[plane] = num_bins;
      parameters_bin_edges[plane] = bin_edges;
    }

    void set_lookup_tables(size_t plane, std::vector<float> tables)
    {
      lookup_tables[plane] = tables;
    }

    size_t digitize(float value, std::vector<float> bin_edges)
    {
      if (value <= bin_edges[0])
        return 0;
      for(size_t i=0; i<bin_edges.size(); i++)
      {
        if (value >= bin_edges[i])
          continue;
        else
          return i-1;
      }
      return bin_edges.size()-1;
    }

    size_t findLookupIndex(float dedx_value, std::vector<float> par_value, size_t plane)
    {
      //findParameterBin
      std::vector<size_t> this_parameters_bins;
      for(size_t i=0; i<par_value.size(); i++)
      {
        size_t aux_index = digitize(par_value[i], parameters_bin_edges[plane][i]);
        this_parameters_bins.push_back(aux_index);
      }

      //findLookUpRow
      size_t lookup_row=0, accumulator_par_bins=1;
      for(size_t i=this_parameters_bins.size(); i-- > 0; )
      {
        lookup_row += (accumulator_par_bins * this_parameters_bins[i]);
        accumulator_par_bins *= parameters_num_bins[plane][i];
      }

      //findLookUpRowindex
      size_t lookup_row_index;
      lookup_row_index = lookup_row * dedx_num_bins[plane];

      //findLookUpRowDedxIndex
      size_t lookup_index = lookup_row_index;
      lookup_index += digitize(dedx_value, dedx_bin_edges[plane]);

      std::cout << "lookup index " << lookup_index << std::endl;
      return lookup_index;
    }

    float LLR_one_hit_one_plane(float dedx_value, std::vector<float> par_value, size_t plane)
    {
      size_t index = findLookupIndex(dedx_value, par_value, plane);
      return lookup_tables[plane][index];
    }

    float LLR_many_hits_one_plane(std::vector<float> dedx_values, std::vector<std::vector<float>> par_values, size_t plane)
    {
      float ll_out = 0;
      for(size_t i=0; i<dedx_values.size(); i++)
      {
        std::vector<float> aux_par;
        for(std::vector<float> par_value: par_values)
        {
          aux_par.push_back(par_value[i]);
        }
        ll_out += LLR_one_hit_one_plane(dedx_values[i], aux_par, plane);
      }
      return ll_out;
    }

  private:
    size_t dedx_num_bins[3];
    std::vector<float> dedx_bin_edges[3];

    std::vector<size_t> parameters_num_bins[3];
    std::vector<std::vector<float>> parameters_bin_edges[3];

    std::vector<float> lookup_tables[3];
  };
}

#endif
