#include "bwlabel.hpp"

BwLabel::BwLabel(){;}
BwLabel::~BwLabel(){;}

//frechet_compare_mat is the input image.
std::vector<std::vector<int>> BwLabel::find_label(cv::Mat input_image)
{
  //Stores the group number of the pixel
  std::vector<std::vector<int>> group_matrix;
  //Matrix containing all the group numbers used
  std::vector<std::vector<int>> used_groups_num;
  
  used_groups_num = first_grouping(input_image, group_matrix );

  second_grouping(used_groups_num, group_matrix);

  final_grouping(used_groups_num, group_matrix);

  return group_matrix;
}


//Goes through each pixels grouping the pixels the first time.
//If the previous pixel is 1 then the group number of previous pixel is taken.
//O group for dark pixels.
std::vector<std::vector<int>> BwLabel::first_grouping(cv::Mat input_image, std::vector<std::vector<int>>& group_matrix)
{
  size_t length_row = input_image.rows;
  size_t length_column = input_image.cols;
  std::vector<int> temp_group_mat;
  std::vector<std::vector<int>> used_groups_num;
  int group_num = 0;
  for (size_t index_row = 0; index_row < length_row; index_row++)
    {
      cv::Scalar intensity = input_image.at<uchar>(index_row,0); 

      int val_zero_col = intensity.val[0];
      if ( val_zero_col== 255) {
	group_num++;
	//if new group_number is assigned it is added to vector group_num_mat
	used_groups_num.push_back({ group_num });
	temp_group_mat.push_back(group_num);
      }
      else {
	temp_group_mat.push_back(0);
      }
      for (size_t index_column = 1; index_column < length_column; index_column++)
	{
	  cv::Scalar cur_intensity =input_image.at<uchar>(index_row,index_column);
	  cv::Scalar prev_intensity = input_image.at<uchar>(index_row,index_column-1);
	  int val_curr_col = cur_intensity.val[0];
	  int val_prev_col = prev_intensity.val[0];
	  //if current pixel is 1 but previous pixel is 0 then increace group number and assign to current pixel
	  if ( val_curr_col == 255 && val_prev_col == 0 )
	    {
	      group_num++;
	      used_groups_num.push_back({ group_num });
	      temp_group_mat.push_back(group_num);
	    }
	  //if previous pixel and current pixel has value 1 use same group number of previous pixel.
	  else if ( val_curr_col == 255 && val_prev_col == 255)
	    {
	      temp_group_mat.push_back(group_num);
	    }
	  //If the current pixel is dark assign group number 0.
	  else if (val_curr_col == 0)
	    {
	      temp_group_mat.push_back(0);
	    }
	}
      
      //each push creates different row in fretchet_group_mat
      group_matrix.push_back(temp_group_mat);
      temp_group_mat.clear();
    }

  //returns matrix containing group number of pixels.
  return used_groups_num;
}


//(all the group numbers used, input vector)
void BwLabel::second_grouping(std::vector<std::vector<int>>& used_groups_num, std::vector<std::vector<int>>& group_matrix )
{
  size_t length_row = group_matrix.size();//number of rows in input img
  size_t length_column = group_matrix[0].size();//number of cols in input img
  size_t length_num = used_groups_num.size();  //Total number of group assigned

  //For each row after first row send the row for processing
  for (size_t index_row = 1; index_row < length_row; index_row++)
    {
      //sends first and second row for processing (ie 0&1, 1&2.........) and vector containing label info.
      process_rows(group_matrix[index_row - 1], group_matrix[index_row], used_groups_num);
    }
}


//Processes previous row and current row sent from group_mat_process function
//(previous row, current row, vector containing all the used groups)
void BwLabel::process_rows(std::vector<int>& first_row, std::vector<int>& sec_row, std::vector<std::vector<int>>& used_groups_num )
{
  
  size_t length_column = first_row.size();
  int last_num = 0;
  //iterate through previous row.
  for (size_t index_column = 0; index_column < length_column; index_column++)
    {
      //if the group number of pixels in previous column is not zero
      if (first_row[index_column]) {
	//vector containing neighbours of bright pixel from previous row in second row.
	std::vector<int> neighbor_vector = find_neighbor(index_column, sec_row);
	for (int neighbor_flag : neighbor_vector)
	  {
	    //making sure that -ve index is not included.
	    if (neighbor_flag != -1) {
	      /*if (index_column >= neighbor_flag)*/
	      {
		//curr_row[neighbour_flag]=group number of neighbour index.
		int neb_group_no = sec_row[neighbor_flag];
		//current group number
		int cur_group_no = first_row[index_column];
		//last member parent of current row
		int cur_las_paren = used_groups_num[sec_row[neighbor_flag] - 1].back();
		//last member parent of previous row
		int prev_las_paren = used_groups_num[first_row[index_column] - 1].back();

		//If there was just one element
		//If the size is 1 or less there is no chance that the last group number in the current row is greater than previous row.
		//If the last group number in the current row is greater than previous row then update the group number to the least one.
		if ((used_groups_num[neb_group_no - 1].size() <= 1) || (cur_las_paren > prev_las_paren))
		  {
		    //stores the last group number of previous index to current
		    used_groups_num[neb_group_no - 1].push_back(used_groups_num[cur_group_no - 1].back());
		  }
		//If the last group number in the current row is less than previous row 
		else if (cur_las_paren< prev_las_paren)
		  {
		    //
		    trace_parent(used_groups_num, used_groups_num[cur_group_no - 1].back(), used_groups_num[neb_group_no - 1].back());
		    used_groups_num[cur_group_no - 1].push_back(used_groups_num[neb_group_no - 1].back());
		    
		  }
	      }
	    }
	  }
      }
    }
}


//finds neighbour of pixel in a row
//(index of current pixel, and current row(not previous) since we are working on previous row)
std::vector<int> BwLabel::find_neighbor(size_t index, std::vector<int>& neighbor_row)
{
  std::vector<int> result;
  //special case if the supplied index is zero.
  if (index == 0)
    {
      if (neighbor_row[0])
	{
	  result.push_back(0);
	}
      if (neighbor_row[1])
	{
	  result.push_back(1);
	}
    }
  //special case if the supplied index is last pixel
  else if (index == neighbor_row.size() - 1)
    {
      if (neighbor_row[neighbor_row.size() - 2])
	{
	  result.push_back(neighbor_row.size() - 2);
	}
      else if (neighbor_row[neighbor_row.size() - 1])
	{
	  result.push_back(neighbor_row.size() - 1);
	}
    }
  else
    {
      if (neighbor_row[index - 1])
	{
	  result.push_back(index - 1);
	}
      else if (neighbor_row[index])
	{
	  result.push_back(index);
	}
      else if (neighbor_row[index + 1])
	{
	  result.push_back(index + 1);
	}
    }
  //returns vector containing connected neighbours.
  return result;
}


//is called when parent group number of current group is less than that of previous
//(vector containing groups used, group number of previous row, group number of current row);
void BwLabel::trace_parent(std::vector<std::vector<int>>& used_groups_num, int pre_num, int tar_num)
{
  //Iterate over all the group number of previous row.
  for (int i = 0; i < used_groups_num[pre_num - 1].size(); i++)
    {
      //for every group number in previous row entry, add new parent. since tar_num is super parent of all.
      int num = used_groups_num[pre_num - 1][i];
      used_groups_num[num - 1].push_back(tar_num);
    }
}


//(vector containing info about group number, matrix containing group number of each pixel)
void BwLabel::final_grouping(std::vector<std::vector<int>>& used_groups_num, std::vector<std::vector<int>>& group_matrix )
{
  size_t length_num = used_groups_num.size();
  size_t length_row = group_matrix.size();
  size_t length_column = group_matrix[0].size();
  //stores curr_num and trace_num
  std::vector<cv::Point2i> traces;
  //iteritate back in vector containing group numbers.
  for (int index_num = length_num - 1; index_num >= 0; index_num--)
    {
      if(used_groups_num[index_num].size() != 1){
      int trace_index = index_num;
      //curr_num = first parent group in current index.
      int curr_num = used_groups_num[trace_index][0], trace_num = curr_num;

      //Iteritate until it reaches to the head of branch
      while (used_groups_num[trace_index].size() != 1)
	{
	  //trace_num = last parenttrace_index
	  trace_num = used_groups_num[trace_index].back();
	  //trace index = index of last parent.
	  trace_index = used_groups_num[trace_index].back() - 1;
	  //sends the vector containing parents of parent vector.
	  remove_duplicate_parent(used_groups_num[trace_index]);
	  //loop ends when it reaches to the top parent.
	}
     
      //will be group number of all connected points.
      trace_num = used_groups_num[trace_index].back();
      //vector containing curr_num and corresponding trace num.
      traces.push_back(cv::Point(curr_num,trace_num));
      }
    }


  for(size_t index_row = 0; index_row <length_row; index_row++)
    {
      for(size_t index_column = 0; index_column <length_column; index_column++)
	{
	  if(group_matrix[index_row][index_column]){
	  for(int i=0; i<traces.size();i++){
	    int curr_num = traces[i].x;
	    int trace_num = traces[i].y;
	    if (group_matrix[index_row][index_column] == curr_num)
	      {
		group_matrix[index_row][index_column] = trace_num;
		break;
	      }	
	  }
	  }
	}
    }
	  
      
  for (int index_row = length_row - 1; index_row > 0; index_row--)
    {
      for (int index_column = length_column - 2; index_column > 0; index_column--)
	{
	  if (group_matrix[index_row][index_column])
	    {

	      std::vector<int> neighbor_vector = find_neighbor(index_column, group_matrix[index_row - 1]);
	      for (int neighbor_flag : neighbor_vector)
		{
		  //if the neighbour from upper row has parent greater than lower row. change to lower one 
		  if (group_matrix[index_row - 1][neighbor_flag]>=group_matrix[index_row][index_column])
		    {
		      group_matrix[index_row - 1][neighbor_flag] = group_matrix[index_row][index_column];
		    }
		  else 
		    {
		      fix_upper(group_matrix, group_matrix[index_row][index_column], group_matrix[index_row - 1][neighbor_flag]);
		    }
		}

	      if (group_matrix[index_row][index_column - 1]>group_matrix[index_row][index_column])
		{
		  group_matrix[index_row][index_column - 1] = group_matrix[index_row][index_column];
		}
	      if (group_matrix[index_row][index_column + 1]>group_matrix[index_row][index_column])
		{
		  group_matrix[index_row][index_column + 1] = group_matrix[index_row][index_column];
		}
	    }
	}
	}
}


//removes the repetation of same parent
void BwLabel::remove_duplicate_parent(std::vector<int>& parents)
{
  std::vector<int> new_array;
  int last_num = parents[0];
  new_array.push_back(last_num);
  for (size_t i = 0; i < parents.size(); i++)
    {
      if (parents[i] != last_num)
	{
	  last_num = parents[i];
	  bool add_flag = true;
	  for (size_t i = 0; i < new_array.size(); i++)
	    {
	      if (last_num == new_array[i])
		{
		  add_flag = false;
		}
	    }
	  if (add_flag)
	    {
	      new_array.push_back(last_num);
	    }
	}
    }
  parents = new_array;
}


//Apply when upper connected neighbour has smaller group number than lower.
//(matrix containing group number, group number from lower column, group number from upper column) 
void BwLabel::fix_upper(std::vector<std::vector<int>>& used_groups_num, int origin_num, int target_num)
{
  size_t length_row = used_groups_num.size();
  size_t length_column = used_groups_num[0].size();
  for (int index_row = 0; index_row < length_row; index_row++)
    {
      for (int index_column = 0; index_column < length_column; index_column++)
	{
	  if (used_groups_num[index_row][index_column] == origin_num)
	    {
	      used_groups_num[index_row][index_column] = target_num;
	    }
	}
    }
}

