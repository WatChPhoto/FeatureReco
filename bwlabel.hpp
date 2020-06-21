#ifndef BWLABEL_HPP
#define BWLABEL_HPP

#include<vector>
#include<opencv2/core.hpp>

class BwLabel{
public:
  BwLabel();
  ~BwLabel();
  std::vector<std::vector<int>> find_label(cv::Mat input_image);
private:
  std::vector<std::vector<int>> first_grouping(cv::Mat input_image, std::vector<std::vector<int>>& group_matrix );
  void second_grouping(std::vector<std::vector<int>>& used_groups_num, std::vector<std::vector<int>>& group_matrix);
  void process_rows(std::vector<int>& first_row, std::vector<int>& sec_row,std::vector<std::vector<int>>& used_groups_num);
  std::vector<int> find_neighbor(size_t index, std::vector<int>& neighbor_row );
  void trace_parent(std::vector<std::vector<int>>& used_groups_num, int pre_num, int tar_num);
  void final_grouping( std::vector<std::vector<int>>& used_groups_num, std::vector<std::vector<int>>& group_matrix );
  void remove_duplicate_parent(std::vector<int>& parents );
  void fix_upper(std::vector<std::vector<int>>& used_groups_num, int orgin_num, int target_num );
 


};

#endif
