%/* ************************************************** */
%     Copyright 2014 Arash Bakhtiari
%
%     you may not use this file except in compliance with the License.
%     You obtain a copy of the License in the LICENSE file
%
%     Unless required by applicable law or agreed to in writing, software
%     distributed under the License is distributed on an "AS IS"" BASIS,
%     WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
%     See the License for the specific language governing permissions and
%     limitations under the License.
%/* ************************************************** */

classdef tree_data < handle
%UNTITLED Summary of this class goes here
%   Detailed explanation goes here

    properties
    end

    methods (Static)
        %/* ************************************************** */
        function [xx,yy,vv] = grid_points(tree)
            xx = []; yy = [];
            cleaves = tree.leaves();
            for lvcnt = 1:length(cleaves)
                cleaf = cleaves{lvcnt};
                if isempty(cleaf.data), continue; end;
                % GRID POINTS
                resPerNode = cleaf.data.resolution;
                [xr,yr,zr,dx,dy,dz] = cleaf.mesh(resPerNode);
                xx = [xx; xr(:)];
                yy = [yy; yr(:)];
            end
        end

        %/* ************************************************** */
        function [vv] = grid_data(tree)
            vv = [];
            cleaves = tree.leaves();
            if isempty(cleaves{1}.data),
                dim=1;
            else
                dim = cleaves{1}.data.dim;
     end;

     for dcnt = 1:dim
         vvtmp = [];
         for lvcnt = 1:length(cleaves)
             cleaf = cleaves{lvcnt};
             if isempty(cleaf.data), continue; end;
             % GRID VALUES
             vals = cleaf.data.values;
             tmp = vals(:,:,:,dcnt);
             vvtmp = [vvtmp; tmp(:)];
         end
         if ~isempty(vvtmp)
             vv(:,dcnt) = vvtmp;
         end
     end
 end

 %/* ************************************************** */
 function init_data(tree, func, resPerNode, t)
     if nargin < 4, t = 0; end;

     cleaves = tree.leaves();
     for lvcnt = 1:length(cleaves)
         cleaf = cleaves{lvcnt};
         [xr,yr,zr,dx,dy,dz] = cleaf.mesh(resPerNode);
         cvalues = func(t,xr,yr,zr);
         cleaf.data.dim          = size(cvalues,4);
         cleaf.data.resolution   = resPerNode;
         for dcnt = 1:size(cvalues,4)
             cleaf.data.values(:,:,:,dcnt)   = cvalues(:,:,:,dcnt);
         end
     end
 end

 %/* ************************************************** */
 function interp(src_tree, dst_tree)
 % interpolate the values of grid points in the destination tree
 % from the values of correspoding nodes in the source tree
 % TODO: this is a brute-force algorithm. optimize it.
     global verbose;

     dst_leaves  = dst_tree.leaves();
     src_leaves  = src_tree.leaves();

     resPerNode  = src_leaves{1}.data.resolution;
     data_dim    = src_leaves{1}.data.dim;

     for dst_lvcnt = 1:length(dst_leaves)
         dst_leaf = dst_leaves{dst_lvcnt};

         % get the mesh of destination leaf
         [dst_xxr,dst_yyr,dst_zzr,dst_dx,dst_dy,dst_dz] = dst_leaf.mesh(resPerNode);

         % init the destination leaf data values
         dst_leaf.data.dim = data_dim;
         dst_leaf.data.resolution = resPerNode;
         % remove the 1 for 3D implementation
         dst_leaf.data.values = zeros([size(dst_xxr) 1 data_dim]);

         for dimcnt = 1:data_dim
             % interpolate the destination values from corresponding source leaves
             tmpval = zeros(size(dst_xxr));
             for src_lvcnt =1:length(src_leaves)
                 src_leaf = src_leaves{src_lvcnt};

                 % find the points inside current source leaf
                 indices = tree_data.points_in_node(src_leaf, dst_xxr, dst_yyr);
                 if ~any(any(indices)), continue; end;

                 xx = dst_xxr(indices);
                 yy = dst_yyr(indices);

                 [src_xxr,src_yyr,src_zzr,src_dx,src_dy,src_dz] = src_leaf.mesh(resPerNode);

                 interp_data = src_leaf.data.values(:,:,:,dimcnt);
                 vv = interp2(src_xxr,src_yyr,interp_data,xx,yy);
                 tmpval(indices) = vv;
             end
             dst_leaf.data.values(:,:,:,dimcnt) = tmpval;
         end
     end
 end

 %/* ************************************************** */
 function val = interp_points(src_tree,xq,yq,zq,INTERP_TYPE)
     src_leaves  = src_tree.leaves();
     resPerNode  = src_leaves{1}.data.resolution;
     data_dim    = src_leaves{1}.data.dim;

     val = zeros([size(xq) 1 data_dim]);
     for dimcnt = 1:data_dim
         % interpolate the destination values from corresponding source leaves
         tmpval = zeros(size(xq));
         for src_lvcnt =1:length(src_leaves)
             src_leaf = src_leaves{src_lvcnt};

             % find the points inside current source leaf
             indices = tree_data.points_in_node(src_leaf, xq, yq);
             if ~any(any(indices)), continue; end;

             xx = xq(indices);
             yy = yq(indices);

             [src_xxr,src_yyr,src_zzr,src_dx,src_dy,src_dz] = src_leaf.mesh(resPerNode);

             interp_data = src_leaf.data.values(:,:,:,dimcnt);
             vv = interp2(src_xxr,src_yyr,interp_data,xx,yy,INTERP_TYPE);
             tmpval(indices) = vv;
         end
         val(:,:,:,dimcnt) = tmpval;
     end
 end

 %/* ************************************************** */
 function max_err = compute_error(tree,fexact, t,INTERP_TYPE)
     leaves  = tree.leaves();
     resPerNode  = leaves{1}.data.resolution;
     data_dim    = leaves{1}.data.dim;
     max_err = 0;
     %for dimcnt = 1:data_dim
     for lvcnt =1:length(leaves)
         leaf = leaves{lvcnt};
         [xr,yr,zr,dx,dy,dz] = leaf.mesh(resPerNode);

         % compute the center of the local grid cells
         xxc = xr+dx/2;
         yyc = yr+dy/2;
         zzc = zr+dz/2;
         xxcc = xxc(1:end-1,1:end-1);
         yycc = yyc(1:end-1,1:end-1);
         zzcc = zzc(1:end-1,1:end-1);

         vale = fexact(t,xxcc,yycc,zzcc);

         valt = tree_data.interp_points(leaf,xxcc,yycc,zzcc,INTERP_TYPE);% leaf.data.values(:,:,:,dimcnt);
                                                                         % compute interpolation error
         diff = vale - valt;
         err = max(max(abs(diff)));
         if err > max_err, max_err = err; end;
     end
     %end
 end

 %/* ************************************************** */
 function [tree_out] = collapse(tree_in)
 % TODO: check that treecells have the exact same structure
 %       -> mids of all treecells are the same.
 % clone structure of the resulting tree from the input trees
     num_trees = length(tree_in);
     num_leaves = length(tree_in{1}.leaves());
     tree_out = qtree.clone(tree_in{1});

     % get the leaves of input trees
     tree_in_leaves = cell(num_trees,num_leaves);
     for tree_in_cnt =1:num_trees
         tree_in_leaves(tree_in_cnt,:) = tree_in{tree_in_cnt}.leaves();
     end

     % iterate over leaves of tree out
     tree_out_leaves = tree_out.leaves();
     for leafcnt = 1:length(tree_out_leaves)
         leaf = tree_out_leaves{leafcnt};

         tmp = tree_in_leaves{1,leafcnt};
         leaf.data.dim           = num_trees;
         leaf.data.resolution    = tmp.data.resolution;
         % TODO: remove the one after extending the code to 3D
         leaf.data.values        = zeros([size(tmp.data.values) 1 num_trees]);

         for tree_in_cnt = 1:num_trees
             tree_in_leaf = tree_in_leaves{tree_in_cnt,leafcnt};
             leaf.data.values(:,:,:,tree_in_cnt) = tree_in_leaf.data.values;
         end
     end
 end

 %/* ************************************************** */
 function plot_grid(tree)
     MS='MarkerSize';
     [txx,tyy] = tree_data.grid_points(tree);
     tree.plottree;
     hold on;
     plot(txx(:),tyy(:),'ro',MS,1);
     axis off; axis equal;
 end

 %/* ************************************************** */
 function plot_data(tree,dim)
     if nargin < 2, dim = 1; end;
     MS='MarkerSize';
     %tree.plottree(0.5)
     hold on
     [txx,tyy]   = tree_data.grid_points(tree);
     [tvv]       = tree_data.grid_data(tree);
     if ~isempty(tvv) & ~isempty(txx) & ~isempty(tyy)
         scatter3(txx,tyy,tvv(:,dim),ones(size(txx)),tvv(:,dim),'filled')
     end
     %plot3(txx(:),tyy(:),tvv(:,dim),'.',MS,1)%,ones(size(txx)),tvv(:,dim),'filled')
     axis off;
     axis equal;
 end

 %/* ************************************************** */
 function indices = points_in_node(node, xx, yy)
 % complication for points that lie right on
 % the boundaries
     [xmin,xmax,ymin,ymax]=corners(node);
     idx = find(xmin <= xx & xx <= xmax);
     idy = find(ymin <= yy & yy <= ymax);
     indices = intersect(idx, idy);
 end
    end
end
