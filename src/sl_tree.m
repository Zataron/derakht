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

function cnext = sl_tree(c,u,v,t,fdo_refine,fconc_exact,fvel_exact)
    global RES_PER_NODE;
    global VERBOSE;
    global INTERP_TYPE;

    VPREVTSTEP  = 1;
    VCURTSTEP   = 2;
    VNEXTSTEP   = 3;
    V2NEXTSTEP  = 4;

    % % VELOCITY (TIME-DEPENDENT)
    % nt = length(t);
    % ucells = cell(1,nt);
    % vcells = cell(1,nt);
    %
    % for tcnt = 1:nt
    %     utmptree = qtree;
    %     utmptree.insert_function(fvelx,fdo_refine,t(tcnt));
    %     qdata.init_data(utmptree,fvelx,RES_PER_NODE,t(tcnt));
    %     ucells{tcnt} = utmptree;
    %
    %     vtmptree = qtree;
    %     vtmptree.insert_function(fvely,fdo_refine,t(tcnt));
    %     qdata.init_data(vtmptree,fvely,RES_PER_NODE,t(tcnt));
    %     vcells{tcnt} = vtmptree;
    % end
    %
    % % MERGE VELOCITY TREES
    % fprintf('-> merge velocity trees\n');
    % num     = size(ucells,2);
    % um_tree = merge(ucells);
    % umcells = clone(um_tree,num);
    % vm_tree = merge(vcells);
    % vmcells = clone(vm_tree,num);
    %
    % % INTERPOLATE VELOCITY VALUES ON THE MERGED TREES
    % fprintf('-> interpolate velocity values on the merged tree\n');
    % for i =1:num
    %     qdata.interp(ucells{i}, umcells{i});
    %     qdata.interp(vcells{i}, vmcells{i});
    % end
    %
    % u = qdata.collapse(umcells);
    % v = qdata.collapse(vmcells);

    %fconc_interp    = fconc_exact;
    %fconc_interp    = @conc_interp;
    fconc_interp     = @conc_interp_qmsl;

    %fvel_interp     = fvel_exact;
    %fvel_interp     = @vel_interp;
    fvel_interp     = @vel_interp_time_independent;

    %fsemilag        = fconc_exact;   
    fsemilag        = @semilag;


    % FIRST METHOD: CONSTRUCT THE NEXT TIME STEP TREE FROM SCRATCH WITH
    %               SEMILAG SOLVER AS REFINEMENT FUNCTION
    % cnext = qtree;
    % cnext.insert_function(fsemilag,fdo_refine);
    % qdata.init_data(cnext,fsemilag,RES_PER_NODE);

    % SECOND METHOD: USE THE PREVIOUS TIME STEP TREE AS STARTING POINT
    %                REFINE/COARSEN WHENEVER IS NEEDED
    % cnext = qtree.clone(c);
    % update_tree(cnext,fdo_refine);

    % THIRD METHOD: SEPARATE THE COMPUTATION AND REFINEMENT PROCESSES
    cnext = qtree.clone(c);

    fprintf('--> SL: ');
    sl_time = tic;
    qdata.init_data(cnext,fsemilag,RES_PER_NODE);
    fprintf('%f\n', toc(sl_time))

    global ADAPTIVE
    if ADAPTIVE
        fprintf('--> Refine: ');
        rt_time = tic;
        refine_tree(cnext, fdo_refine);
        fprintf('%f\n',toc(rt_time))
    end

    %/* ************************************************** */
    function do_coarsen = refine_tree(node, fvisit)
        do_coarsen = true;
        kidsval = true;
        if ~node.isleaf
            for k=1:4
                kidval = refine_tree(node.kids{k}, fvisit);
                kidsval = kidsval & kidval;
            end
            if kidsval
                [refine_node] = fvisit(node,fsemilag,t(VNEXTSTEP));
                if refine_node, do_coarsen = false; return;
                else
                    % COARSEN NODE
                    if VERBOSE,
                        fprintf('--> coarsen node: ')
                        fprintf(' level %2d: anchor:[%1.4f %1.4f] \n', ...
                                node.level,node.anchor(1),node.anchor(2));
                    end
                    coarsen(node)
                    qdata.set_node_fn(node, fsemilag, RES_PER_NODE, t);
                    do_coarsen = true;
                    % do_coarsen = false;
                    return;
                end
            else do_coarsen = false; return;
            end
        else
            [refine_node] = fvisit(node,fsemilag,t(VNEXTSTEP));
            if refine_node
                % REFINE NODE
                if VERBOSE,
                    fprintf('--> refine node: ')
                    fprintf(' level %2d: anchor:[%1.4f %1.4f] \n', ...
                            node.level,node.anchor(1),node.anchor(2));
                end
                refine(node);
                for kcnt=1:4,
                    qdata.set_node_fn(node.kids{kcnt}, fsemilag, RES_PER_NODE, t);
                end;
                for kcnt=1:4, refine_tree(node.kids{kcnt},fvisit); end;
                do_coarsen = false;
                return;
            else
                % set_node_values(node);
                do_coarsen = true;
                return;
            end
        end

    end

    % %/* ************************************************** */
    % function val = update_tree(node, fvisit)
    %     val = true;
    %     kidsval = true;
    %     if ~node.isleaf
    %         for k=1:4
    %             kidval = update_tree(node.kids{k}, fvisit);
    %             kidsval = kidsval & kidval;
    %         end
    %         if kidsval
    %             [refine_node, values] = fvisit(node,fsemilag,t(VNEXTSTEP));
    %             if refine_node
    %                 val = false;
    %                 return;
    %             else
    %                 % COARSEN THE NODE
    %                 if VERBOSE,
    %                     fprintf('--> coarsen node: ')
    %                     fprintf(' level %2d: anchor:[%1.4f %1.4f] \n', ...
    %                             node.level,node.anchor(1),node.anchor(2));
    %                 end
    %                 coarsen(node)
    %                 set_node_values(node, values);
    %                 val = true;
    %                 return;
    %             end
    %         else
    %             val = false;
    %             return;
    %         end
    %     else
    %         [refine_node, values] = fvisit(node,fsemilag,t(VNEXTSTEP));
    %         if refine_node
    %             % REFINE THE NODE
    %             if VERBOSE,
    %                 fprintf('--> refine node: ')
    %                 fprintf(' level %2d: anchor:[%1.4f %1.4f] \n', ...
    %                         node.level,node.anchor(1),node.anchor(2));
    %             end
    %             refine(node);
    %             for kcnt=1:4, update_tree(node.kids{kcnt},fvisit); end;
    %             val = false;
    %             return;
    %         else
    %             set_node_values(node, values);
    %             val = true;
    %             return;
    %         end
    %     end

    %     %/* ************************************************** */
    %     function set_node_values(node, values)
    %         if VERBOSE,
    %             fprintf('--> set semilag values for node: ')
    %             fprintf(' level %2d: anchor:[%1.4f %1.4f] \n', ...
    %                     node.level,node.anchor(1),node.anchor(2));
    %         end
    %         node.data.dim           = 1;
    %         node.data.resolution    = RES_PER_NODE;
    %         node.data.values        = values;
    %     end
    % end

    %/* ************************************************** */
    function refine(node)
        node.data = [];
        node.create_kids();
    end

    %/* ************************************************** */
    function coarsen(node)
        node.kids = [];
        node.isleaf = true;
    end

    %/* ************************************************** */
    %function val = QMSL(tdummy,x,y,z)
    %    if VERBOSE, fprintf('--> calling QMSL function!\n'); end
    %    val = semilag_rk2(x,y,z,fconc_interp,fvel_interp,t);
    %end

    %/* ************************************************** */
    %function valq = quasi_monotone_adjust(tdummy,x,y,z,valq)
    %    valq = qdata.QMSL_adjust(c,x,y,z,valq,INTERP_TYPE);        
    %end    

    %/* ************************************************** */
    function val = semilag(tdummy,x,y,z)
        if VERBOSE, fprintf('--> calling semilag function!\n'); end
        val = semilag_rk2(x,y,z,fconc_interp,fvel_interp,t);
    end

    %/* ************************************************** */
    function ci = conc_interp(tq,xq,yq,zq)
        ci = qdata.interp_points(c,xq,yq,zq,INTERP_TYPE);
        ci = conc_out(ci,xq,yq,zq);
    end

    %/* ************************************************** */
    function ci = conc_interp_qmsl(tq,xq,yq,zq)
        ci = qdata.interp_points_qmsl(c,xq,yq,zq,INTERP_TYPE);
        ci = conc_out(ci,xq,yq,zq);
    end

    %/* ************************************************** */
    function cq = conc_out(cq,xq,yq,zq)
        % OUTSIDE THE SIMULATION DOMAIN
        %ce = fconc_exact(tq,xq,yq,zq);
        ce = zeros(size(cq));
        out = xq<0 | xq>1  | yq<0 | yq>1 | zq<0 | zq>1;
        cq(out) = ce(out);
    end

    %/* ************************************************** */
    function [uq,vq,wq] = vel_interp(tq,xq,yq,zq)
        uval = qdata.interp_points(u,xq,yq,zq);
        vval = qdata.interp_points(v,xq,yq,zq);
        [uq,vq,wq,] = interp_vel_temporal(uval,vval,0,t,tq,INTERP_TYPE);
        [uq,vq,wq] = vel_out(uq,vq,wq,xq,yq,zq);

        function [uq,vq,wq] = vel_out(uq,vq,wq,xq,yq,zq)
        % OUTSIDE THE SIMULATION DOMAIN
            out = xq<0 | xq>1  | yq<0 | yq>1 | zq<0 | zq>1;
            %             [ue, ve, we] = fvel_exact(tq,xq,vq,zq);
            %             uq(out) = ue(out);
            %             vq(out) = ve(out);
            %             wq(out) = we(out);
            ue = zeros(size(uq));
            uq(out) = ue(out);
            vq(out) = ue(out);
            wq(out) = ue(out);
        end
    end

    %/* ************************************************** */
    function [uq,vq,wq] = vel_interp_time_independent(tq,xq,yq,zq)
        uval = qdata.interp_points(u,xq,yq,zq,INTERP_TYPE);
        vval = qdata.interp_points(v,xq,yq,zq,INTERP_TYPE);
        uq = uval;
        vq = vval;
        wq = zeros(size(uq));
        [uq,vq,wq] = vel_out(uq,vq,wq,xq,yq,zq);

        function [uq,vq,wq] = vel_out(uq,vq,wq,xq,yq,zq)
        % OUTSIDE THE SIMULATION DOMAIN
            out = xq<0 | xq>1  | yq<0 | yq>1 | zq<0 | zq>1;
            [ue, ve, we] = fvel_exact(tq,xq,vq,zq);
            uq(out) = ue(out);
            vq(out) = ve(out);
            wq(out) = we(out);
        end
    end

    %/* ************************************************** */
    function [mt] = merge(treecells)
        mt = treecells{1};
        for counter=2:length(treecells)
            mt = qtree.merge(mt, treecells{counter});
        end
    end

    %/* ************************************************** */
    function [tree_clones] = clone(tree_src, num)
        tree_clones = cell(1,num);
        for counter=1:num
            tree_clones{counter} = qtree.clone(tree_src);
        end
    end

    function plotvel()
        figure('Name','SEMI-LAG QUAD-TREES');
        subplot(3,4,2)
        c.plottree;
        qdata.plot_data(c);
        title('c(t)');

        subplot(3,4,3)
        cnext.plottree;
        qdata.plot_data(cnext);
        title('c(t+dt)');

        subplot(3,4,5)
        ucells{1}.plottree;
        qdata.plot_data(ucells{1});
        title('u(t(n-1))');

        subplot(3,4,6)
        ucells{2}.plottree;
        qdata.plot_data(ucells{2});
        title('u(t(n))');

        subplot(3,4,7)
        ucells{3}.plottree;
        qdata.plot_data(ucells{3});
        title('u(t(n+1))');

        subplot(3,4,8)
        ucells{4}.plottree;
        qdata.plot_data(ucells{4});
        title('u(t(n+2))');

        subplot(3,4,9)
        vcells{1}.plottree;
        qdata.plot_data(vcells{1});
        title('v(t(n-1))');

        subplot(3,4,10)
        vcells{2}.plottree;
        qdata.plot_data(vcells{2});
        title('v(t(n))');

        subplot(3,4,11)
        vcells{3}.plottree;
        qdata.plot_data(vcells{3});
        title('v(t(n+1))');

        subplot(3,4,12)
        vcells{4}.plottree;
        qdata.plot_data(vcells{4});
        title('v(t(n+2))');
    end
end
