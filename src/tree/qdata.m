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

classdef qdata < handle
%UNTITLED Summary of this class goes here
%   Detailed explanation goes here

    properties
    end

    methods (Static)

        %/* ************************************************** */
        function init_data(tree, func, resPerNode, t)
        % Initializes the values of all leaf octants based on the
        % input function values and grid resolution per node
            if nargin < 4, t = 0; end;
            cleaves = tree.leaves();
            for lvcnt = 1:length(cleaves)
                cleaf = cleaves{lvcnt};
                qdata.set_node_fn(cleaf, func, resPerNode, t);
            end
        end              
        
        %/* ************************************************** */
        function init_data_alt(tree, func, resPerNode, t)
        % Initializes the values of all leaf octants based on the
        % input function values and grid resolution per node
        % computes all at once rather than going leaf by leaf
        % NOTE: only works for regular grid
            if nargin < 4, t = 0; end;
            
            % gather all grid points in [xx,yy]
            [xx,yy] = qdata.grid_points_alt(tree);
            zz = zeros(size(xx));  
            
            % invoke func  
            vv = func(t,xx,yy,zz);
            
            % set values for all leaves
            pos = 1;
            points_per_node = (resPerNode+1)*(resPerNode+1);
            cleaves = tree.leaves();             
            for lvcnt = 1:length(cleaves)
                cleaf = cleaves{lvcnt};
                v = reshape(vv(pos:pos+points_per_node-1, 1),[resPerNode+1,resPerNode+1]);
                qdata.set_node_val(cleaf, v, resPerNode);
                pos = pos + points_per_node;
            end
        end          

        %/* ************************************************** */
        function interp_tree(src_tree, dst_tree)
        % interpolate the values of grid points in the destination tree
        % from the values of correspoding nodes in the source tree
            global VERBOSE;
            global INTERP_TYPE;
            global RES_PER_NODE;
            dst_leaves  = dst_tree.leaves();
            for dst_lvcnt = 1:length(dst_leaves)
                dst_leaf = dst_leaves{dst_lvcnt};
                [xx,yy,zz,dx,dy,dz] = dst_leaf.mesh(RES_PER_NODE, INTERP_TYPE);
                tmpval = qdata.interp_points(src_tree, xx, yy, zz, INTERP_TYPE);
                qdata.set_node_val(dst_leaf, tmpval, RES_PER_NODE);
            end
        end
        
        %/* ************************************************** */
        function vv = QMSL_adjust(src_leaf,xx,yy,interp_data,vv,dx,dy)
            % QMSL stage
            % xx,yy are departure points in src_leaf
            % interp_data is src_leaf.data.values
            % vv is the high order solution
            [xmin,xmax,ymin,ymax] = corners(src_leaf);  
            for i = 1:length(xx)
                % get values at the corners of the box around current x,y
                x = xx(i);
                y = yy(i);
                if (x == xmax || y == ymax)
                    continue;
                end;
                x_index = 1+floor((x-xmin)/dx);
                y_index = 1+floor((y-ymin)/dy);

                M = interp_data(y_index:y_index+1, x_index:x_index+1);

                % get min, max and adjust
                upper = max(M(:));
                lower = min(M(:));
                val = vv(i);
                if val < lower
                    vv(i) = lower;
                elseif val > upper
                    vv(i) = upper;
                end
            end                       
        end
        
        %/* ************************************************** */
        function valq = interp_points(src_tree,xq,yq,zq,INTERP_TYPE)
            global RES_PER_NODE;
            global ENABLE_QMSL;
            src_leaves  = src_tree.leaves();
            valq = zeros(size(xq));
            for src_lvcnt =1:length(src_leaves)
                src_leaf = src_leaves{src_lvcnt};
                indices = qdata.points_in_node(src_leaf, xq, yq);
                if ~any(any(indices)), continue; end;
                xx = xq(indices);
                yy = yq(indices);
                if strcmp(INTERP_TYPE, 'CHEBYSHEV')
                    global CHEB_IMPL
                    if strcmp(CHEB_IMPL, 'IAS')
                        vv = zeros(size(xx));
                        w = src_leaf.data.values;
                        [xmin,xmax,ymin,ymax] = src_leaf.corners;
                        xs = (xx - xmin)*2/(xmax-xmin)-1.0;
                        ys = (yy - ymin)*2/(ymax-ymin)-1.0;
                        for xindx =1:size(xx,1)
                            vv(xindx) =  cheb.chebeval2(w, xs(xindx), ys(xindx));
                        end
                    elseif strcmp(CHEB_IMPL, 'CHEBFUN')
                        w = src_leaf.data.values;
                        vv = w(xx,yy);
                    end
                else
                    % Regular Grid
                    [xxr,yyr,zzr,dx,dy,dz] = src_leaf.mesh(RES_PER_NODE);
                    interp_data = src_leaf.data.values;
                    vv = interp2(xxr,yyr,interp_data,xx,yy, INTERP_TYPE);                   
                    
                    % QMSL adjustment
                    if ENABLE_QMSL
                        vv = qdata.QMSL_adjust(src_leaf,xx,yy,interp_data,vv,dx,dy);
                    end
                end
                valq(indices) = vv;
            end
        end            
        
        %/* ************************************************** */
        function valq = interp_points_cqmsl(src_tree,xq,yq,zq,S,INTERP_TYPE)
            %NOTE: xq,yq must be column vectors
            global RES_PER_NODE;
            global T0_MASS;
            valq = zeros(size(xq));
            if strcmp(INTERP_TYPE, 'CHEBYSHEV')
                % CQMSL only supported for regular grid
                return;
            end
            valq_low = valq;
            valq_high = valq;             
            src_leaves = src_tree.leaves();
            for src_lvcnt =1:length(src_leaves)                 
                src_leaf = src_leaves{src_lvcnt};
                indices = qdata.points_in_node(src_leaf, xq, yq);
                if ~any(any(indices)), continue; end;
                xx = xq(indices);
                yy = yq(indices);
                
                % High order and low order solution
                [xxr,yyr,zzr,dx,dy,dz] = src_leaf.mesh(RES_PER_NODE);
                interp_data = src_leaf.data.values;
                vv = interp2(xxr,yyr,interp_data,xx,yy, INTERP_TYPE);
                valq_high(indices) = vv;
                valq_low(indices) = interp2(xxr,yyr,interp_data,xx,yy, 'linear');

                % QMSL solution
                valq(indices) = qdata.QMSL_adjust(src_leaf,xx,yy,interp_data,vv,dx,dy);
            end
            
            % Conservation stage
            %calculate mass difference               
            mass_difference = sum(sum(S.*valq))-T0_MASS;
            if mass_difference == 0
                return;
            end 
            %calculate weights
            weight = max(0,sign(mass_difference)*((valq_high-valq_low).^3));
            if ~any(any(weight))
                return;
            end
            %calculate lambda
            lambda = mass_difference/sum(sum(weight.*S));
            %calculate new valq
            valq = valq - lambda*weight;          
        end        

        %/* ************************************************** */
        function max_err = compute_error(tree, fexact, t,INTERP_TYPE)
            leaves     = tree.leaves();
            resPerNode = leaves{1}.data.resolution;
            data_dim   = leaves{1}.data.dim;
            max_err = 0;
            max_err_pnt = [];
            max_err_node = [];
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
                valt = qdata.interp_points(leaf,xxcc,yycc,zzcc,INTERP_TYPE);

                diff = vale - valt;
                [err, indx] = max(abs(diff(:)));
                [i_row, i_col] = ind2sub(size(diff),indx);
                if err > max_err
                    vale;
                    valt;
                    max_err = err;
                    max_err_node = leaf;
                    max_err_pnt = [i_row, i_col];
                end;
            end
        end

        % %/* ************************************************** */
        % function [tree_out] = collapse(tree_in)
        % % TODO: check that treecells have the exact same structure
        % %       -> mids of all treecells are the same.
        % % clone structure of the resulting tree from the input trees
        %     num_trees = length(tree_in);
        %     num_leaves = length(tree_in{1}.leaves());
        %     tree_out = qtree.clone(tree_in{1});

        %     % get the leaves of input trees
        %     tree_in_leaves = cell(num_trees,num_leaves);
        %     for tree_in_cnt =1:num_trees
        %         tree_in_leaves(tree_in_cnt,:) = tree_in{tree_in_cnt}.leaves();
        %     end

        %     % iterate over leaves of tree out
        %     tree_out_leaves = tree_out.leaves();
        %     for leafcnt = 1:length(tree_out_leaves)
        %         leaf = tree_out_leaves{leafcnt};

        %         tmp = tree_in_leaves{1,leafcnt};
        %         leaf.data.dim           = num_trees;
        %         leaf.data.resolution    = tmp.data.resolution;
        %         % TODO: remove the one after extending the code to 3D
        %         leaf.data.values        = zeros([size(tmp.data.values) 1 num_trees]);

        %         for tree_in_cnt = 1:num_trees
        %             tree_in_leaf = tree_in_leaves{tree_in_cnt,leafcnt};
        %             leaf.data.values(:,:,:,tree_in_cnt) = tree_in_leaf.data.values;
        %         end
        %     end
        % end

        %/* ************************************************** */
        function set_node_fn(node, fn, resPerNode, t)
            global INTERP_TYPE
            [xx,yy,zz,dx,dy,dz] = node.mesh(resPerNode, INTERP_TYPE);
            fnval = fn(t,xx,yy,zz);
            qdata.set_node_val(node, fnval, resPerNode);
        end

        %/* ************************************************** */
        function set_node_val(node, fnval, resPerNode)
            global INTERP_TYPE
            if strcmp(INTERP_TYPE, 'CHEBYSHEV')
                w = qdata.get_node_cheb_interpolant(node, fnval, resPerNode);
                node.data.values = w;
            else
                node.data.values = fnval;
            end
            node.data.dim = 1;
            node.data.resolution = resPerNode;
        end

        %/* ************************************************** */
        function [w] = get_node_cheb_interpolant(node, fn_val, resPerNode)
            global CHEB_IMPL;
            if strcmp(CHEB_IMPL, 'IAS')
                w = cheb.chebcoeff(fn_val);
            elseif strcmp(CHEB_IMPL, 'CHEBFUN')
                [xmin xmax ymin ymax] = node.corners;
                w = chebfun2(fn_val, [xmin xmax ymin ymax]);
            end
        end

        %/* ************************************************** */
        function [val] = get_node_values_cheb(node, resPerNode)
            global CHEB_IMPL;
            if strcmp(CHEB_IMPL, 'IAS')
                x = cheb.chebnodes1(resPerNode);
                y = x;
                w = node.data.values;
                val = cheb.chebeval2(w,x,y);
            elseif strcmp(CHEB_IMPL, 'CHEBFUN')
                w = node.data.values;
                [xx,yy,zz,dx,dy,dz] = node.mesh(resPerNode, ...
                                                'CHEBYSHEV');
                val = w(xx,yy);
            end
        end

        %/* ************************************************** */
        function [xx,yy,vv] = grid_points(tree)
            global INTERP_TYPE;
            global RES_PER_NODE;
            xx = []; yy = [];
            cleaves = tree.leaves();
            for lvcnt = 1:length(cleaves)
                cleaf = cleaves{lvcnt};
                if isempty(cleaf.data), continue; end;
                [xr,yr,zr,dx,dy,dz] = cleaf.mesh(RES_PER_NODE,INTERP_TYPE);
                xx = [xx; xr(:)];
                yy = [yy; yr(:)];
            end
        end
        
        %/* ************************************************** */
        function [xx,yy,vv] = grid_points_alt(tree)
            % difference to grid_points(): also returns grid points of leaves without data
            global INTERP_TYPE;
            global RES_PER_NODE;
            xx = []; yy = [];
            cleaves = tree.leaves();
            for lvcnt = 1:length(cleaves)
                cleaf = cleaves{lvcnt};
                %if isempty(cleaf.data), continue; end;
                [xr,yr,zr,dx,dy,dz] = cleaf.mesh(RES_PER_NODE,INTERP_TYPE);
                xx = [xx; xr(:)];
                yy = [yy; yr(:)];
            end
        end         

        %/* ************************************************** */
        function [vv] = grid_data(tree)
            global INTERP_TYPE;
            vv = [];
            cleaves = tree.leaves();
            for lvcnt = 1:length(cleaves)
                cleaf = cleaves{lvcnt};
                if isempty(cleaf.data), continue; end;
                % GRID VALUES
                if strcmp(INTERP_TYPE,'CHEBYSHEV')
                    %evaluate chebyshev polynomials at the cheb
                    resPerNode = cleaf.data.resolution;
                    vals = qdata.get_node_values_cheb(cleaf, resPerNode);
                    tmp = vals(:,:);
                else
                    vals = cleaf.data.values;
                    tmp = vals(:,:,:);
                end
                vv = [vv; tmp(:)];
            end
        end

        %/* ************************************************** */
        function plot_grid(tree)
            MS='MarkerSize';
            [txx,tyy] = qdata.grid_points(tree);
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
            [txx,tyy]   = qdata.grid_points(tree);
            [tvv]       = qdata.grid_data(tree);
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
        
        %/* ************************************************** */
        function val = get_mass(src_tree,INTERP_TYPE)
            global RES_PER_NODE;
            if strcmp(INTERP_TYPE, 'CHEBYSHEV') 
                % CHEBYSHEV GRID
                val = 0;
                src_leaves  = src_tree.leaves();
                for src_lvcnt =1:length(src_leaves)
                    src_leaf = src_leaves{src_lvcnt};     
                    
                    % calculate length*width for each quadrant
                    [xx,yy,zz,dx,dy,dz] = src_leaf.mesh(RES_PER_NODE,INTERP_TYPE);
                    xs = zeros(size(xx));
                    ys = zeros(size(yy));
                    for i = 1:(size(xx,2) - 1)
                        xs(:,i) = xx(:,i+1) - xx(:,i);
                    end
                    for i = 1:(size(yy,1) - 1)
                        ys(i,:) = yy(i+1,:) - yy(i,:);
                    end     
                    H = ys .* xs;
                    
                    % get value for the anchor of each quadrant
                    vs = qdata.get_node_values_cheb(src_leaf, RES_PER_NODE);
                          
                    % add up length*width*value for each quadrant
                    val = val + sum(sum(vs .* H));
                end
                return;
            end
            
            % REGULAR GRID
            val = 0;
            src_leaves  = src_tree.leaves();
            for src_lvcnt =1:length(src_leaves)
                src_leaf = src_leaves{src_lvcnt};
                interp_data = src_leaf.data.values;
                [xxr,yyr,zzr,dx,dy,dz] = src_leaf.mesh(RES_PER_NODE);
                
                x = xxr(1,1:end);
                y = yyr(1:end,1);
                val = val + trapz(y,trapz(x,interp_data,2));
            end
        end                            
        
        %/* ************************************************** */
        function val = get_mass_squared(src_tree,INTERP_TYPE)
            if strcmp(INTERP_TYPE, 'CHEBYSHEV')
                % currently only regular grid is supported!
                val = 0;
                return;
            end
            global RES_PER_NODE;
            src_leaves  = src_tree.leaves();
            val = 0;
            for src_lvcnt =1:length(src_leaves)
                src_leaf = src_leaves{src_lvcnt};
                interp_data = src_leaf.data.values;
                [xxr,yyr,zzr,dx,dy,dz] = src_leaf.mesh(RES_PER_NODE);
                
                x = xxr(1,1:end);
                y = yyr(1:end,1);
                val = val + trapz(y,trapz(x,interp_data.^2,2));
            end
        end   
        
        %/* ************************************************** */
        function val = get_height(src_tree,INTERP_TYPE)
            % get highest point in simulation domain
            % only used for testing purposes
            if strcmp(INTERP_TYPE, 'CHEBYSHEV')
                % currently only regular grid is supported!
                val = 0;
                return;
            end
            src_leaves  = src_tree.leaves();
            val = -1;
            for src_lvcnt =1:length(src_leaves)
                interp_data = src_leaves{src_lvcnt}.data.values;
                interp_max = max(interp_data(:));
                if (val < interp_max)
                    val = interp_max;
                end
            end
        end         
        
        %/* ************************************************** */
        function [e_diss, e_disp, e_sum, e_total] = get_interpolation_errors(src_tree, fexact, t)
            global INTERP_TYPE;
            
            [X,Y] = qdata.grid_points(src_tree);
            [interp_values] = qdata.grid_data(src_tree);    
            real_values = fexact(t,X,Y,0);

            %get covariance, standard deviations, means, correlation coeff
            cv = cov(real_values, interp_values);
            std_real = std(real_values);
            std_interp = std(interp_values);
            mean_real = mean(real_values);
            mean_interp = mean(interp_values);
            corr = cv(1,2)/(std_real*std_interp);
            
            %calculate errors
            e_diss = (std_real - std_interp)^2 + (mean_real - mean_interp)^2;
            e_disp = 2*(1-corr)*std_real*std_interp;  
            e_sum = e_diss+e_disp;
            
            %e_total is the mean squared error
            %it should be equal to e_sum
            M = (real_values-interp_values).^2;
            e_total = sum(M(:))/length(M(:));
        end    

        %/* ************************************************** */
        function draw_cross_section(src_tree,cs,mode,fig_num)
            if nargin < 4, fig_num = 1; end;
            global INTERP_TYPE;
            if strcmp(INTERP_TYPE, 'CHEBYSHEV')                    
                if strcmp(mode,'cross_X')
                    yy = reshape(linspace(0,1,99),[],1);
                    xx = cs*ones(size(yy));    
                    ii = yy;
                else
                    xx = reshape(linspace(0,1,99),[],1);
                    yy = cs*ones(size(xx));  
                    ii = xx;
                end                
                  
                zz = zeros(size(xx));
                iv = qdata.interp_points(src_tree,xx,yy,zz,INTERP_TYPE);
                
                figure(fig_num)
                plot(ii, iv)
                grid on   
            else
                if strcmp(mode,'cross_X')
                   [X,Y] = qdata.grid_points(src_tree);
                   [V] = qdata.grid_data(src_tree);
                   [iy,I] = sort(Y(X == cs));

                   iv = V(X == cs);
                   iv = iv(I);

                   figure(1)
                   plot(iy, iv)
                   grid on  
                else  
                   [X,Y] = qdata.grid_points(src_tree);
                   [V] = qdata.grid_data(src_tree);
                   [ix,I] = sort(X(Y == cs));

                   iv = V(Y == cs);
                   iv = iv(I);

                   figure(1)
                   plot(ix, iv)
                   grid on               
                end
            end
        end     
        
    end
end
