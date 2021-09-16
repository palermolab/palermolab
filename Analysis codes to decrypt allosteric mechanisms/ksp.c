#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <float.h>
#include <limits.h>
#include <string.h>
#define number_of_res 1362
#define k_paths_wanted 10
float dist[number_of_res][number_of_res];
int path[number_of_res];
float d[number_of_res];
int continous_path[number_of_res];
int Alist[k_paths_wanted][number_of_res];
int continous_path_pointer=0;
int org_continous_path[number_of_res];
int temp_path[number_of_res];
float temp_distance=INT_MAX;
int candidate_paths[number_of_res][number_of_res];
float candidate_paths_distance[number_of_res];
  int candidate_path_number=-1;
  float root;
/* The algorithm being implemented here can be written as follows:
 * A[1] := Shortest path from S to T
for k := 2 to K do
    for i := 1 to length(A[1])-1 do
        nodeA = A[k-1].node(i)
        for j := 1 to k-2
            nodeB = A[j].node(i)
            if (nodeA == nodeB) then
                distance(nodeA,nodeB) = infinity
            endif
        endfor
        S[i] = The shortest path from nodeA to T according to the current distance values
        R[i] = The path in A[k-1] from S to nodeA
        B[i] = R[i] + S[i]
    endfor
    A[k] = min length paths amongst all B[i]
    restore distance(nodeA,nodeB) to original value if modified
endfor
*/
void stack_push(int residue_to_push)
{
       continous_path[++continous_path_pointer]=residue_to_push;
}

void dijiktra(int source,int destination)
{
int i,j;
int mini;
int counter=0;
int visited[number_of_res];
for(i=0;i < number_of_res;i++)
{
 for(j=0;j < number_of_res;j++)
 {
   d[i]=INT_MAX;
   visited[i]=0;
   path[i]=number_of_res;
 }
}
d[source]=0;

for(i=0;i < number_of_res;i++)
{
  mini=-1;
 for(j=0;j < number_of_res;j++)
	if((visited[j]!= 1) && ((mini==-1) || d[j] < d[mini]))
		mini=j;
 
 visited[mini]=1;
 if(mini != destination)
 {
 for(j=0;j < number_of_res;j++)
 {
   if(dist[mini][j]!=0)
	    if(d[mini]+dist[mini][j] < d[j])
	    {
		d[j]= d[mini]+dist[mini][j];
		path[j]=mini;
	    }
  }
}
}
  
}

int distance_finder(float spur_distance, int current_path[number_of_res], int intermediate)
{
  float distance_to_check=0;
  float precision = 0.00001;
  int i=0;
  root=0;
//    printf("Mutating %d\n\n",intermediate);
  while(current_path[i]!=intermediate)
  {
    distance_to_check=distance_to_check+dist[current_path[i]][current_path[i-1]];
//     printf(" Adding up %d %d Value is %f\n",current_path[i], current_path[i-1],distance_to_check);
    ++i;
  }
  distance_to_check=distance_to_check+dist[current_path[i]][current_path[i-1]];
//   printf("Adding up %d %d Value is %f\n",current_path[i], current_path[i-1],distance_to_check);
  root=distance_to_check;
  distance_to_check=distance_to_check+spur_distance;
//   printf("After adding spur %f\n\n\n",distance_to_check);
  for(i=0;i <= (candidate_path_number-1);++i)
  {
//     printf("Comparing %f %f\n",candidate_paths_distance[i],distance_to_check);
    if(((candidate_paths_distance[i]-precision)< distance_to_check) && ((candidate_paths_distance[i] + precision)> distance_to_check))
    {
//       printf("YES\n");
      return 1;
    }
  }  
  return 0;
}

void add_to_candidate(int candidate_path_number,int k_th_path, int destination,float temp_distance_val)
{
  int i=0;
  while(continous_path[i]!= destination)
  {
    candidate_paths[candidate_path_number][i]=continous_path[i];
//     printf("%d ",candidate_paths[candidate_path_number][i]);
    ++i;
  }
  candidate_paths[candidate_path_number][i]=continous_path[i];
//    printf("%d\n",candidate_paths[candidate_path_number][i]);
}


int path_comparison(int start, int source, int destination, int k, int current_path[number_of_res])
{
  int i=source;
  int j=0;
  int flag=0;
  int tempnode;
  int new_path=0;
  float temp_value_distance;
  float temp[number_of_res][number_of_res];
  for(j=0; j<=k-2;++j) //comparing all previous paths with k-1th path which are in Alist[0] thru Alist[k-2]
  {
    flag=0;
    for(i=0; i <=(source);i++)//comparing upto the ith node,called the source, of the jth path with the current_path
    {

	  if(Alist[j][i]==current_path[i])
	  {
	    flag=flag+1;
	  }
    }
	  if(flag == source + 1)//this flag is only true if the first i nodes in both cases match
	  {
	    dist[Alist[j][source]][Alist[j][source+1]]=dist[Alist[j][source+1]][Alist[j][source]]=0;
	  }
  }
	    dist[current_path[source]][current_path[source+1]]=dist[current_path[source+1]][current_path[source]]=0;
//  	    printf("Calling on %d %d\n",current_path[source],destination);
	    dijiktra(current_path[source],destination);
	    if((d[destination]!=INT_MAX) && (distance_finder(d[destination],current_path,current_path[source])==0))
	    {
	      	    ++candidate_path_number;
		    concatinate(start,current_path[source],source,destination,current_path);
		    	    
	    candidate_paths_distance[candidate_path_number]=root+d[destination];
// 	    printf("Number of Candidate Paths is %d\n",candidate_path_number);
	    add_to_candidate(candidate_path_number,k,destination,temp_distance);
	    }
}

int remove_top_path(int destination)
{
  int i=0;
  int j=0;
  if(candidate_path_number==0)
  {
    --candidate_path_number;
  }
  else
  {    
  for(i=1;i <=candidate_path_number;i++)
  {
      candidate_paths_distance[i-1]=candidate_paths_distance[i];
      j=0;
      while(candidate_paths[i][j]!=destination)
      {
	candidate_paths[i-1][j]=candidate_paths[i][j];
	++j;
      }
      candidate_paths[i-1][j]=candidate_paths[i][j];
  }
  --candidate_path_number;
  }
  return 0;
}



swap_paths(int path_number,int destination)
{
  int i=0;
  int temp_path[number_of_res];
//   printf("\n\nSwapping Starts Here\n\n");
  while(candidate_paths[path_number][i]!=destination)
  {
    temp_path[i]=candidate_paths[path_number][i];
//     printf("%d ",temp_path[i]);
    ++i;
  }
  temp_path[i]=candidate_paths[path_number][i];
//       printf("%d\n",temp_path[i]);
  
  i=0;
  while(candidate_paths[path_number+1][i]!=destination)
  {
    candidate_paths[path_number][i]=candidate_paths[path_number+1][i];
//     printf("%d ",candidate_paths[path_number][i]);
    ++i;
  }
  candidate_paths[path_number][i]=candidate_paths[path_number+1][i];
//   printf("%d \n",candidate_paths[path_number][i]);
   
    i=0;
  while(temp_path[i]!=destination)
  {
    candidate_paths[path_number+1][i]=temp_path[i];
//     printf("%d ",candidate_paths[path_number+1][i]);
    ++i;
  }
    candidate_paths[path_number+1][i]=temp_path[i];
//     printf("%d\n",candidate_paths[path_number+1][i]);
    
//  printf("Swapping Ends Here\n\n");
}



sort_candidate_paths(int candidate_path_number,int destination)
{
  int i;
  int bubble_sort_flag=1;
  float temp_swap;
  while(bubble_sort_flag==1)
  {
  bubble_sort_flag=0;
  for(i=0; i < (candidate_path_number);++i)
    if(candidate_paths_distance[i] > candidate_paths_distance[i+1])
    {
      bubble_sort_flag=1;
      swap_paths(i,destination);
      temp_swap=candidate_paths_distance[i];
      candidate_paths_distance[i]=candidate_paths_distance[i+1];
      candidate_paths_distance[i+1]=temp_swap;
    }
  }
}

       
concatinate(int source, int intermediate,int current_node, int destination,int current_path[number_of_res])
{
  continous_path_pointer=-1;
  int i=0;
//   printf("Problem is here %d\n",intermediate);
  while(current_path[i]!=intermediate)
  {
  stack_push(current_path[i]);
  ++i;
  }
  stack_push(current_path[i]);
  pathfinder(intermediate,destination);
  stack_push(destination);
}
       
ksp(float dist[number_of_res][number_of_res],int source,int destination,int k_paths,int current_node, int max_pointer,int current_path[number_of_res])
{
  float temp[number_of_res][number_of_res];
  float temp_val;
  float op_d=0;
  int flag=0;
  int i=0;
  int j=0;
  for(i=0;i<number_of_res;++i)//making a copy of the matrix
  {
    for(j=0;j<number_of_res;++j)
    {
      temp[i][j]=dist[i][j];
    }
  }

   for(current_node=0;current_node <=(max_pointer-1);++current_node)//for every node
   {
    path_comparison(source,current_node,destination,k_paths,current_path);//check is
   	for(i=0;i<number_of_res;++i)//restoring the matrix
	      {
	for(j=0;j<number_of_res;++j)
		{
		  dist[i][j]=temp[i][j];
		}
	      }
   }
	  sort_candidate_paths(candidate_path_number,destination);
// 	    i=0;
// 		printf("\n\nPRINTING ALL SORTED PATHS\n\n");
// 	  for(j=0;j <= (candidate_path_number) ;++j)
// 	  {
// 	  i=0;
// 	  while(candidate_paths[j][i]!=destination)
// 	  {
// 	    printf("%d ",candidate_paths[j][i]);
// 	    ++i;
// 	  }
// 	  printf("%d --%f\n",candidate_paths[j][i],candidate_paths_distance[j]);
// 	  }
// 	  printf("\n\nEND OF ALL SORTED PATHS\n\n");

	  i=0;
// 	  printf("The %dnd shortest path is:",k_paths);
	  while(candidate_paths[0][i]!=destination)
	  {
	    continous_path[i]=candidate_paths[0][i];
	    printf("%d ",continous_path[i]);
	    ++i;
	  }
	  continous_path[i]=candidate_paths[0][i];
	  temp_distance=candidate_paths_distance[0];
	  printf("%d	%f\n",continous_path[i],temp_distance);
	  remove_top_path(destination);
	  
	  i=0;
// 	 printf("\n\nAfter Removal PRINTING ALL SORTED PATHS\n\n");
// 	  for(j=0;j <=(candidate_path_number) ;++j)
// 	  {
// 	  i=0;
// 	  while(candidate_paths[j][i]!=destination)
// 	  {
// 	    printf("%d ",candidate_paths[j][i]);
// 	    ++i;
// 	  }
// 	  printf("%d --%f\n",candidate_paths[j][i],candidate_paths_distance[j]);
// 	  }
// 	  printf("\n\nEND OF ALL SORTED PATHS AFTER REMOVAL\n\n");

	  
	  
      return 0;
}




reverse_pathfinder(int source,int destination)
{
  if(path[destination]==number_of_res)
     {
      return 0;
     }
     else if (path[destination]==source)//direct path has been found so we don't need to call the function again. 
     {
     }
     else
     {
     stack_push(path[destination]);
     //printf("%d ",path[destination]);
     reverse_pathfinder(source,path[destination]);//getting the intermediate nodes starting from the first one
     //push the node to the top i.e. continous_path[0]==source and continous_path[N]=destination
    }
return continous_path_pointer;
}


pathfinder(int source,int destination)
{
  if(path[destination]==number_of_res)
     {
      return 0;
     }
     else if (path[destination]==source)//direct path has been found so we don't need to call the function again. 
     {
     }
     else
     {
     pathfinder(source,path[destination]);//getting the intermediate nodes starting from the first one
     stack_push(path[destination]);//push the node to the top i.e. continous_path[0]==source and continous_path[N]=destination
     //printf("%d ",path[destination]);
    }
return continous_path_pointer;
}

    
push_dijiktra(int source, int destination)
{
  int temp;
  int temp_node;
  continous_path_pointer=-1;

    // printf("Candidate :%d ",source);
     stack_push(source);//pushing the source onto the stack using stack_push
     pathfinder(source,destination);
     stack_push(destination);//pushing the destination onto the stack using stack_push
    //printf("%d --%f\n",destination,d[destination]);
}



int main()
{
  int a, i, j, n;
  char line[8192];
  char * pch;
  {
FILE *fin = fopen("dist.txt", "r");
        fscanf(fin, "%d", &n);
        for (i = 0; i < number_of_res; ++i)
                for (j = 0; j < number_of_res; ++j)
                {
                        fscanf(fin, "%f", &dist[i][j]);
                }
        fclose(fin);

for (i=0; i< number_of_res; ++i) //this loop calculates every log
{
  for(j=0; j< number_of_res; ++j)
      dist[i][j]=fabsl(log10(dist[i][j]));
}

FILE *four = fopen("notfour.txt", "r");
  i=0;
    if(four==NULL)
	 {
	  printf("ERROR:The notfour.txt file does not exist.\n");
	  return 0;
	}
    else
	{
	  while( i < number_of_res)//only the first number_of_res lines will be read, others will be disregared
	  {
	  if(fgets(line, 8192, four)!=NULL)//gets a new line for each residue/index i
	  {
	    pch =strtok(line,"\n");//since the newline character is included in fgets we need to get rid of it. the pch has the entire line with spaces except for the new line
	    pch = strtok (pch," ");//get the first number from the string list stored in the pch
	    while (pch != NULL)//if the number is not empty(end of line)
		{
		    a=atoi(pch);//convert to an integer
		    a=a-1;//subtract 1 to get 0 indexing
		if ((a != i+1) && ( a != i-1))//if the number is not the next neighbor of index i
		{
		    dist[i][a]=0;//setting the distance to 0
		    dist[a][i]=0;//setting the distance to 0
		}
		else 
		{
		  
		}
		    pch = strtok (NULL," ");//get the next number
		}
	  }
	      i++;//get the next index
	  }
	}
fclose(four);
  }
int temp[number_of_res];
int k;
int current_path[number_of_res];
int start_node=X;
int destination_node=X;
// for(start_node=0;start_node<number_of_res;++start_node)
// {
//   for(destination_node=(start_node+1);destination_node<number_of_res;++destination_node)
//   {
    candidate_path_number=-1;
    dijiktra(start_node,destination_node);
    push_dijiktra(start_node,destination_node);
    i=0;
    while (continous_path[i]!=destination_node)
    {
      printf("%d ",continous_path[i]);
      ++i;
    }
    printf("%d	%f\n",continous_path[i],d[destination_node]);
    for(k=2;k<=k_paths_wanted;k++)
    {
      i=0;
    while(continous_path[i]!=destination_node)
    {
      Alist[k-2][i]=continous_path[i];
      current_path[i]=continous_path[i];
      ++i;
    }
    Alist[k-2][i]=continous_path[i];
    current_path[i]=continous_path[i];
    ksp(dist,start_node,destination_node,k,0,i,current_path);
    }
//   }
// }
}
  


