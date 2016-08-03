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

classdef cheb < handle
%UNTITLED Summary of this class goes here
%   Detailed explanation goes here

    properties
    end

    methods (Static)

        %/* ************************************************** */
        function [x] = chebnodes1(N, xmin, xmax)
            if nargin < 2, xmin=-1; xmax=1; end;
            i = [0:N]';
            global CHEB_KIND
            if CHEB_KIND == 1
                x = -cos((i+1/2)*pi/(N+1));
            elseif CHEB_KIND == 2
                x = -cos(i*pi/N);
            end
            if nargin < 2, return; end;
            x = 0.5*(xmin+xmax) + 0.5*(xmax-xmin).*x;
        end

        %/* ************************************************** */
        function [xx, yy] = chebnodes2(resX, resY, xmin, xmax, ymin, ...
                                        ymax)
            if nargin < 2, resY = resX; end;
            if nargin < 3, xmin=-1; xmax=1; ymin=-1; ymax=1;end;
            x = cheb.chebnodes1(resX, xmin, xmax);
            y = cheb.chebnodes1(resY, ymin, ymax);
            [xx, yy] = meshgrid(x,y);
        end

        function [w] = chebcoeff(fn_val)
            global CHEB_KIND
            n1 = size(fn_val, 2);
            n2 = size(fn_val, 1);
            x = cheb.chebnodes1(n1-1);
            y = cheb.chebnodes1(n2-1);
            Tx = cheb.chebpoly(n1-1,x);
            Ty = cheb.chebpoly(n2-1,y);

            if CHEB_KIND == 1
                Nx = n1; 
                Ny = n2;
            elseif CHEB_KIND == 2
                Nx = n1-1;
                Ny = n2-1;
                fn_val(1,:) = fn_val(1,:)/2;
                fn_val(:,1) = fn_val(:,1)/2;
                fn_val(end,:) = fn_val(end,:)/2;
                fn_val(:,end) = fn_val(:,end)/2;
            end

            % using discrete orthogonality of chebyshev polynomials to compute the
            % coefficients for approximation using chebyshev
            % polynomial basis.
            for i=1:n1
                w_(:,i) = (reshape(fn_val(:,i),1,[])*Ty)*2/Ny;
            end
            for j=1:n2
                w(j,:) = (reshape(w_(j,:),1,[])*Tx)*2/Nx;
            end
        end
        
        %/* ************************************************** */
        function fval = filterfun(ord,x)
            if ord == 0
                %0 order: Fejer
                fval = 1-x;
            elseif ord == 1
                %first order: Lanczos filter
                fval = sin(pi*x)./(pi*x);
                fval(x == 0) = 1;       
            elseif ord == 2
                %second order: raised cosine
                fval = 0.5*(1+cos(pi*x));
            else
                %construct a filter of order ord
                a = -36.0437;
                fval = exp(-a*x.^ord);
            end
        end
        
        %/* ************************************************** */
        function filter = chebfilter_shifted(N)
            % indices from -N/2 to N/2
            bot = -floor((N-1)/2);
            top = ceil((N-1)/2);           
            
            filter = cheb.filterfun(1, linspace(bot,top,N)/top);
        end
        
        %/* ************************************************** */
        function filter = chebfilter(N)
            % indices from 0 to N-1
            filter = cheb.filterfun(1, linspace(0,N-1,N)/(N-1));
        end        
        
        function coeff = chebfilter2d(n,m,N)
            tau = 0.03;
            r = 0.3;
            %coeff = 1 ./ ( 1 + tau*((n+m).^r).^2 );
            
            coeff = cheb.filterfun(1,(n+m)/(N-2));
        end

        %/* ************************************************** */
        function [fval] = chebeval2(w,x,y)
        % CHEBEVAL2(W, X, Y) Compute the values of a chebyshev
        % approximation at a regular grid specified by X, Y. where
        % W is the corresponding chebyshev coefficients.
            global CHEB_KIND
            global FILTER
            fval = zeros(length(y),length(x));
            n1 = size(w,2);
            n2 = size(w,1);
            T_x = cheb.chebpoly(n1-1, x);
            T_y = cheb.chebpoly(n2-1, y);

            w(1,:) = w(1,:)/2;
            w(:,1) = w(:,1)/2;
            if CHEB_KIND == 2
                w(end,:) = w(end,:)/2;
                w(:,end) = w(:,end)/2;
            end
            
            w_ = w;

            % 1-d filter in y-direction
            if FILTER
             filter_y = cheb.chebfilter( n2 );
             T_y = T_y .* repmat(filter_y, size(T_y,1), 1);
            end
            
            % 1-d filter in x-direction
            if FILTER
             filter_x = cheb.chebfilter( n1 );
             T_x = T_x .* repmat(filter_x, size(T_x,1), 1);
            end
            
            % 2-d filter
%             M_mat = repmat( linspace(0,n1-1,n1), n2, 1 );               % column indices
%             N_mat = repmat( reshape(linspace(0,n2-1,n2),[],1), 1, n1 ); % row indices
%             filter_coeffs = cheb.chebfilter2d(N_mat, M_mat, n1+n2);
%             w_ = w_ .* filter_coeffs;
            
            %NOTES:
            %T_y: length(y) x n2
            %T_x: length(x) x n1
            %f_(:,i): length(y) x 1
            %f_: length(y) x n1
            %fval: length(y) x length(x)
            for i=1:size(w_,2)
                f_(:,i)=T_y*reshape(w_(:,i),[],1);
            end
            for j=1:length(y)
                fval(j,:)=T_x*reshape(f_(j,:),[],1);
            end
        end

        %/* ************************************************** */
        function [T] = chebpoly(n,x)
        % CHEBPOLY(N, X) Retruns values of all chebyshev polynomials upto degree N
        % at points X
            x_ = reshape(x,[],1);
            T0 = ones(size(x_));
            T(:,1) = T0;
            if n == 0
                return;
            end
            T1 = x_;
            T(:,2) = T1;
            if n == 1
                return;
            end
            for i = 2:n
                T(:,i+1) = 2*x_.*T1-T0;
                T0 = T1;
                T1 = T(:,i+1);
            end
        end
    end
end
