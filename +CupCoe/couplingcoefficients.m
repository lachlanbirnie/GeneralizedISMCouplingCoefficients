function [A] = couplingcoefficients(k, Room, Source, Receiver)
% [A] = couplingcoefficients(k, Room, Source, Receiver)
% =========================================================================
% File:         +CupCoe\couplingcoefficients.m
% Creation:     20-09-18
% Author:       Lachlan Birnie
% Version:      7.30 (29-11-19)
%
% Description:
%
%   Solves and returns the generalized image source method coupling 
%   coefficients [alpha_vu_nm(k)], for a given room, source region,
%   receiver region, and a single frequency k.
%
% Inputs:
%
%   k       (double) Wave number (frequency)
%            - units: m^-1
%            - format: scalar (@TODO: or 1 dimensional vector)
%
%   Room    (struct) Properties of rectangular room:
%
%           Room.D      : [1] (integer) depth order of image source method.
%
%           Room.b      : [1 by 6] (double) vector of wall reflection 
%                         coefficients, [bx1 bx2 by1 by2 bz1 bz2]
%
%           Room.size   : [1 by 3] (double) vector of Room dimensions,
%                         [x y z] w.r.t bottom left corner in metres.
%
%   Source  (struct) Properties of source region:
%
%           Source.Os   : [1 by 3] (double) position of the source regions
%                         origin w.r.t room's bottom left corner,
%                         [Xos Yos Zos] in metres.
%
%           Source.Rs   : [1] (double) radial size of the source region,
%                         (Rs) in metres.
%
%           Source.N    : [1] (integer) truncation order of source region,
%                         @NOTE: If not defined, function uses ceil(k*Rs).
%
%   Reciever (struct) Properties of receiver region:
%
%           Receiver.Or : [1 by 3] (double) position of the receiver region
%                         origin w.r.t the room's bottom left corner, 
%                         [Xor Yor Zor] in metres. 
%
%           Receiver.Rr : [1] (double) radial size of the receiver region, 
%                         (Rr) in metres. 
%
%           Reciever.V  : [1] (integer) truncation order of the receiver 
%                         region.
%                         @NOTE: If not defined function uses ceil(k*Rr).
%
% Outputs:
%
%   A       Generalized image source coupling coefficients
%            - [(V+1)^2 by (N+1)^2 {by length(k)}] matrix of coefficients.
%            - index of v^th order u^th mode in 1st dimension is v^2+v+u+1
%            - index of n^th order u^th mode in 2nd dimension is n^2+n+m+1
%              @TODO: k^th dimension still needs to be added. 
%
% Dependencies:
%
%   cupcoe_sphericalharmonicmatrix(q, Harm, Image)
%   cupcoe_sphericalhankelmatrix(q, k, Image)
%   cupcoe_wignermatrix(q, Harm, Image)
%   cupcoe_wigner3jvector(j1, j2, j3, m1, m2, m3)
%
% History Log:
%
%   20-09-18    - function and its i/o is defined.
%   26-09-18    - improved Wigner 3j to solve comlumn instead of matrix.
%   27-09-18    - improved Ynm to have faster mapping to zero padded
%                 matrix.
%               - v4: Ynm is only copied for valid orders now.
%   02-10-18    - Y matrix now uses 'Fourier Acoustics' definition. 
%                 As a result, new coefficients now match the old ones.
%               - Function now checks, and accounts for a concentric case.
%   31-10-18    - Updates to naming conventions.
%               - Project is now in a +CupCoe MATLAB package.
%
%   12-06-19 v7.10  - trying to get point to region coupling coefficients.
%
%   30-09-19 v7.20  - Fixed wigner function valid input criteria. 
%                     The mistake in the selection criteria from previous 
%                     versions did not influence the coupling coefficient 
%                     calculation.
%
%   29-11-19 v7.30  - Fixed math error when sources are positioned on the 
%                     z-axis. Credit goes to Huanyu for finding this issue.
%
% =========================================================================

  %% ---------------------------------------------------------------------
  %     Input Checking
  %  ---------------------------------------------------------------------
    
    % k:
    validateattributes( k ...
                        ,{'double'} ...
                        ,{'nonempty','positive','nonzero','vector'} ...
                        );
    
    % Room:
    if ~isstruct(Room)
        error('Expected input Room to be a struct');
    else
        validateattributes( Room.D ...
                            ,{'double'} ...
                            ,{'scalar','integer'} ...
                            );
        validateattributes( Room.b ...
                            ,{'double'} ...
                            ,{'size',[1 6],'>=',0,'<=',1} ...
                            );
    end
    
    % Source:
    if ~isstruct(Source)
        error('Expected input Source to be a struct');
    else    
        validateattributes( Source.Os ...
                            ,{'double'} ...
                            ,{'2d','vector','size',[1 3],'positive'} ...
                            );    
        validateattributes( Source.Rs ...
                            ,{'double'} ...
                            ...,{'scalar','positive'} ...
                            ,{'scalar'} ... % changed 12-06-19 v7.10
                            );
        if isfield(Source,'N')
            validateattributes( Source.N ...
                                ,{'double'} ...
                                ,{'scalar','>=',0,'integer'} ...
                                );
        else
            Source.N = ceil(max(k(:)) * Source.Rs);
        end
    end
    
    % Receiver: 
    if ~isstruct(Receiver)
        error('Expected input Receiver to be a struct');
    else    
        validateattributes( Receiver.Or ...
                            ,{'double'} ...
                            ,{'2d','vector','size',[1 3],'positive'} ...
                            );    
        validateattributes( Receiver.Rr ...
                            ,{'double'} ...
                            ,{'scalar','positive'} ...
                            );      
        if isfield(Receiver,'V')
            validateattributes( Receiver.V ...
                                ,{'double'} ...
                                ,{'scalar','>=',0,'integer'} ...
                                );
        else
            Receiver.V = ceil(max(k(:)) * Receiver.Rr);
        end
    end
    
    %% Check For Concentric Case - - - - - - - - - - - - - - - - - - - - - 
        
        if (Source.Os == Receiver.Or)
            f_concentric = 1;
            fprintf('couplingcoefficients(): direct path is removed\n');
        else
            f_concentric = 0;
        end
    
  %% ---------------------------------------------------------------------
  %     Image Indexing (Image)
  %  ---------------------------------------------------------------------
  
    % depth indexing (d) [(2D+1)^3 by 1]:
    d1 = repmat((-Room.D:Room.D).', [1, 2*Room.D+1, 2*Room.D+1]);
    d1 = d1(:);
    d2 = repmat((-Room.D:Room.D), [2*Room.D+1, 1, 2*Room.D+1]);
    d2 = d2(:);
    d3 = repmat(permute((-Room.D:Room.D), [3,1,2]), ...
                [2*Room.D+1, 2*Room.D+1, 1]);
    d3 = d3(:);
    
    % wall indexing (p) [1 by 8]:    
    p1 = [0, 1, 0, 1, 0, 1, 0, 1];    
    p2 = [0, 0, 1, 1, 0, 0, 1, 1];
    p3 = [0, 0, 0, 0, 1, 1, 1, 1];
    
    % image indexing [1 by 8(2D+1)^3]:
    %{
    %     - These fields contain the values of the index terms d1,d2,d3
    %       p1,p2,p3 for all of the images, 1 to 8*(2D+1)^3
    % 
    %     Ex: The index terms for the 404^th image are given as:
    %         d(404) = (Image.d1(404),Image.d2(404),Image.d3(404))
    %         p(404) = (Image.p1(404),Image.p2(404),Image.p3(404))
    %
    %
    %     @note: Used length instead of 8 or (2D+1)^3 so that the d and p
    %            vectors above can be changed to image subsets. 
    %}
    Image.d1 = repmat(d1, [1 length(p1)]);
    Image.d1 = Image.d1(:).';
    Image.d2 = repmat(d2, [1 length(p1)]);
    Image.d2 = Image.d2(:).';
    Image.d3 = repmat(d3, [1 length(p1)]);
    Image.d3 = Image.d3(:).';
    Image.p1 = repmat(p1, [length(d1) 1]);
    Image.p1 = Image.p1(:).';
    Image.p2 = repmat(p2, [length(d1) 1]);
    Image.p2 = Image.p2(:).';
    Image.p3 = repmat(p3, [length(d1) 1]);
    Image.p3 = Image.p3(:).';
    
    % check for concentric case:
    if f_concentric
        % remove direct path (middle index, where d = 0 & p = 0):
        dpIndex = ceil((2*Room.D+1)^3 /2);
        Image.d1(dpIndex) = [];
        Image.d2(dpIndex) = [];
        Image.d3(dpIndex) = [];
        Image.p1(dpIndex) = [];
        Image.p2(dpIndex) = [];
        Image.p3(dpIndex) = [];
    end
    
    % count number of images:
    Image.num = length(Image.d1);        
    
    % clear temp variables:
    clear d1 d2 d3 p1 p2 p3 dpIndex; 
    
  %% ---------------------------------------------------------------------
  %     Image Properties (Image)
  %  ---------------------------------------------------------------------
  
    % solve relfection coefficient for each image [1 by 8(2D+1)^3]:
    %{
    %     Ex: The reflection coefficient of the 404^th image is given by:
    %         Image.B(404);
    %}
    Image.B =   (  Room.b(1).^abs(Image.d1 - Image.p1) ...
                .* Room.b(2).^abs(Image.d1)  ...
                .* Room.b(3).^abs(Image.d2 - Image.p2) ...
                .* Room.b(4).^abs(Image.d2) ...
                .* Room.b(5).^abs(Image.d3 - Image.p3) ...
                .* Room.b(6).^abs(Image.d3) ...
                );
            
    % solve image depth vectors (Rd) [8(2D+1)^3 by 3]:
    Image.Rd = [2 .* Image.d1.' .* Room.size(1), ...
                2 .* Image.d2.' .* Room.size(2), ...
                2 .* Image.d3.' .* Room.size(3)];
    
    % solve image wall vectors (Rp) [8(2D+1)^3 by 3]:
    xs = Source.Os(1);
    ys = Source.Os(2);
    zs = Source.Os(3);
    xr = Receiver.Or(1);
    yr = Receiver.Or(2);
    zr = Receiver.Or(3);
    
    Image.Rp = [xr - xs + (2 .* Image.p1.' .* xs), ...
                yr - ys + (2 .* Image.p2.' .* ys), ...
                zr - zs + (2 .* Image.p3.' .* zs)];
            
    % solve image to receiver vectors (Rp+Rd) [1 by 8(2D+1)^3]:
    Image.x = (Image.Rp(:,1) + Image.Rd(:,1)).';
    Image.y = (Image.Rp(:,2) + Image.Rd(:,2)).';
    Image.z = (Image.Rp(:,3) + Image.Rd(:,3)).';
    Image.r = sqrt( Image.x.^2 + Image.y.^2 + Image.z.^2 );
    Image.theta = (pi/2) - atan2(Image.z, sqrt(Image.x.^2 + Image.y.^2));
    Image.phi = wrapTo2Pi( atan2( Image.y, Image.x ));
    
    % clear temp variables:
    clear xs ys zs xr yr zr;
    
  %% ---------------------------------------------------------------------
  %     Harmonic Indexing (Harm)
  %  ---------------------------------------------------------------------
    
    % constants:
    V = Receiver.V;
    N = Source.N;    
    
    % create [(V+1)^2 by (N+1)^2] matrix of harmonic index terms vu, nm:
    v = repmat( repelem((0:V), 2.*(0:V)+1).', [1, (N+1)^2]);
    u = repmat((1:(V+1)^2).', [1, (N+1)^2]) - v.^2 - v - 1;
    n = repmat( repelem((0:N), 2.*(0:N)+1), [(V+1)^2, 1]);
    m = repmat((1:(N+1)^2), [(V+1)^2, 1]) - n.^2 - n - 1;
    
    % save them as [(V+1)^2 * (N+1)^2 by 1] vectors of harmonic terms:
    %{
    %     - These fields provide the harmonic index terms of the i^th 
    %       harmonic coupling coefficient.
    % 
    %     Ex: The harmonic terms of the 5^th coupling coefficient is given:
    %         (v,u) = (Harm.v(5), Harm.u(5))
    %         (n,m) = (Harm.n(5), Harm.m(5))
    %}
    Harm.v = v(:);
    Harm.u = u(:);
    Harm.n = n(:);
    Harm.m = m(:);
    
    % count number of coupled harmonics:
    Harm.num = length(Harm.n);
    
    % clear temp variables:
    clear V N v u n m;

  %% ---------------------------------------------------------------------
  %     Solve Coupling Coefficiens
  %  ---------------------------------------------------------------------
  
    % constants:
    N = Source.N;
    V = Receiver.V;
  
    % define coupling coefficient matrix:
    A = zeros(Harm.num, Image.num);
        
    %% loop over Hankel addition theorm - - - - - - - - - - - - - - - - - 
    %{
    %   - This loop relates to the (L=0 -> n+v) summation in the Snvmu
    %     equation (Eq.14). 
    %
    %   - For this function, the L summation is denoted by 'q', for 
    %     'q' = 0 to N+V.
    % 
    %   - For each loop iteration, the q dependent terms are solved for 
    %     all [(V+1)^2 (N+1)^2 by 8(2D+1)^3] harmonic and image pairs,
    %     and are then added to the running sum of coupling coefficients.
    %}
    
    for q = (0 : N+V)
        
        % Spherical harmonic terms, Y_l_(m' - u) (t, p): 
        qthY = conj( cupcoe_sphericalharmonicmatrix(q, Harm, Image) );
        
        % Spherical Hankel terms, h_l (k Image.r):
        qthH = cupcoe_sphericalhankelmatrix(q, k, Image);
        
        % Wigner terms W1W2Z:
        qthW = cupcoe_wignermatrix(q, Harm, Image);

        % Add q^th term of Hankel addition theorem equation:
        A = A + ((1i)^q .* qthY .* qthH .* qthW);

        % feedback progress:
        fprintf('couplingcoefficients(): q progress = %i / %i\n', q, N+V);
        
    end %q
            
    %% add non-q dependent coefficient terms - - - - - - - - - - - - - - -
    
    %1) (-1)^( (p2+p3)m + p3n) ):
        term1 = (-1).^((Image.p2 + Image.p3).*Harm.m + Image.p3.*Harm.n);
    
    %2) 4*pi*i^(v-n):
        term2 = (4*pi).*(1i).^(Harm.v - Harm.n);
        
    %3) (-1)^(2m' - u), m' = (-1^(p1+p2) *m)
        term3 = (-1).^(2.*((-1).^(Image.p1 + Image.p2).*Harm.m) - Harm.u);
        
    %A) multiply the terms along with reflection coefficients (B):
        A = Image.B .* term1 .* term2 .* term3 .* A;
        
    % @NOTE: A is now [(V+1)^2 (N+1)^2 by 8(2D+1)^3]
        
    %% process final coupling coefficients - - - - - - - - - - - - - - - -
    
        % sum over image sources -> [(V+1)^2 (N+1)^2 by 1]:
        A = sum(A,2);
        
        % reformat -> [(V+1)^2 by (N+1)^2]:
        A = reshape(A, [(Receiver.V+1)^2, (Source.N+1)^2]);
            
end
% end    [A] = couplingcoefficients(k, Room, Source, Receiver)

% ======================================================================= %
%                               FUNCTIONS                                 %
% ======================================================================= %

function [Y] = cupcoe_sphericalharmonicmatrix(q, Harm, Image)
% [Y] = cupcoe_sphericalharmonicmatrix(q, Harm, Image)
% -------------------------------------------------------------------------
% File:     couplingcoefficients.m
% Date:     21-09-18
% Author:   Lachlan Birnie
%
% Description:
%
%   Returns Y_(q)(u-m') (theta, phi) in matrix form, for a given order q 
%   (Hankel addition theorem) and the matrix elements defined by input 
%   Harm and Image strutures
%
%   - m' denotes (-1)^((p1+p2)*m)
%
%   - Y is [(V+1)^2 (N+1)^2 by 8(2D+1)^3]
%
%   Spherical harmonic definition used:
%
%   Y = sqrt( (2n+1)/4pi (n-m)!/(n+m)! ) Pnm(cos t) e^(i m phi) 
%
%   where: Pnm = (-1)^m (1-X^2)^(M/2) * (d/dX)^M { P(N,X) }
%
%   for postive m, and: Pn(-m) = (-1)^m (n-m)!/(n+m)! Pnm
%
% Inputs:
%
%   q       Harmonic order of the complete spherical harmonic matrix [Y]
%
%   Harm    Structure with fields .n .m .v .u, which define the first
%           dimension of the spherical harmonic matrix [Y].
%
%   Image   Structure with fields .p1 .p2 for m', and .theta .phi for 
%           (theta, phi) terms of the spherical harmonic matrix [Y].
%
% Outputs:
%
%   Y       Matrix of spherical harmonic terms of an order 'q' when the 
%           coupling mode indexes |(u - m')| are <= q
%
%           [Y(h,i)] = Y_[ q, Harm.u(h)-m' ] (Image.theta(i), Image.phi(i))
%
%           m' = (-1)^(Image.p1(i) + Image.p2(i)) * Harm.m(h)
%
%   @note:  Because the harmonic order/mode pairs of the input matrix can 
%           be invalid for a Yq(u-m') spherical harmonic term, 
%           we set the elements of Y where |u-m'| > q to zero.
%
%           For efficiency, we also keep the Yq(u-m') terms as zero when:
%           n+v < q, or |n-v| > q, as these terms are zero for the Wigner
%           3jm symbols of the Hankel addition theorem. 
%
%           This allows the matrices Y to be summed together over q,
%           without the need of mapping elements to higher orders.
%
% Change Log:
%
%   26-09-18    - Fixed typo, w = (u - m') [was: (m'-u)]
%               - Fixed indexing issue, matrix is now of order L (q)
%                 and rows of valid modes (u - m') <= q are non-zero. 
%   02-10-18    - Changed spherical harmonic definition to that of the 
%                 'Fourier Acoustics' book. 
%   27-11-18    - Found issue: when q = 86 and g = 85, qthY will have NaN
%                 elements. This is caused from 
%                 sqrt(factorial(q - (-g)) / factorial(q + (-g))) = inf
%                 ( factorial(q - g) ./ factorial(q + g) ) = 0
%                 and, int*0 = NaN. Result should be zero.
%               - HARD CODED FIX @TODO find a proper solution.
%
%% ------------------------------------------------------------------------
    
    % non-valid q(m'-u) indexes remain zero:
    Y = zeros(length(Harm.n), length(Image.theta));
    
    % find columns for [Y] where m' = m, (positive images):
    i_pimg = (-1).^(Image.p1 + Image.p2) == 1;
    % find columns for [Y] where m' = -m, (negative images):
    i_nimg = (-1).^(Image.p1 + Image.p2) == -1;
    
    % solve legendre for order q and all modes 0 -> q (and -q -> -1):
    pvwVector = legendre(q,cos(Image.theta));
    
    for g = 0:q

        % solve Y_(q,g)
        gthY = ( sqrt( (2*q+1) / (4*pi) ) ...
                .* sqrt(factorial(q - g) / factorial(q + g)) ...
                .* pvwVector(abs(g)+1,:) ...
                .* exp(1i .* g .* Image.phi) ...
                );

        % find valid rows of [Y] where (u-m) = g, (positive images):
        i_pmode = ((Harm.u - Harm.m) == g) ...
                & ((Harm.n + Harm.v) >= q) ...
                & (abs(Harm.n-Harm.v) <= q);                
        % find valid rows of [Y] where (u--m) = g, (negative images):
        i_nmode = ((Harm.u -- Harm.m) == g) ...
                & ((Harm.n + Harm.v) >= q) ...
                & (abs(Harm.n-Harm.v) <= q);                       

        % Save Y(q,g) to valid harmonic + image pairs:
        Y(i_pmode, i_pimg) = repmat(gthY(i_pimg), [sum(i_pmode), 1]);
        Y(i_nmode, i_nimg) = repmat(gthY(i_nimg), [sum(i_nmode), 1]);
            
        % Repeat for negative mode (-g):
        if g ~= 0
            
            % TEMPORARY HARD CODE FIX TO HIGH ORDERS BEING NAN (27-11-18)
            if ((factorial(q - g)/factorial(q + g)) ~= 0)
                % solve Y_(q,-g): 
                ngthY = ( sqrt( (2*q+1) / (4*pi) ) ...
                        .* sqrt(factorial(q - (-g)) / factorial(q + (-g))) ...
                        .* (-1)^(g) ...
                        .* ( factorial(q - g) ./ factorial(q + g) ) ...
                        .* pvwVector(abs(g)+1,:) ...
                        .* exp(1i .* -g .* Image.phi) ...
                        );
            else
                ngthY = zeros(size(gthY));
            end
            % END - TEMP FIX
                
            % find valid rows of [Y] where (u-m) == -g:
            i_pmode = ( ((Harm.u - Harm.m) == -g) ...
                    & ((Harm.n + Harm.v) >= q) ...
                    & (abs(Harm.n-Harm.v) <= q) ...
                    );
            % find valid rows of [Y] where (u--m) == -g:
            i_nmode = ( ((Harm.u -- Harm.m) == -g) ...
                    & ((Harm.n + Harm.v) >= q) ...
                    & (abs(Harm.n-Harm.v) <= q) ...
                    );

            % Save Y+(q,0) to valid harmonic,images:
            Y(i_pmode, i_pimg) = repmat(ngthY(i_pimg), [sum(i_pmode), 1]);
            Y(i_nmode, i_nimg) = repmat(ngthY(i_nimg), [sum(i_nmode), 1]);
            
        end %-gth mode
        
    end %gth mode
    
end
% end    [Y] = cupcoe_sphericalharmonicmatrix(q, Harm, Image)

function [H] = cupcoe_sphericalhankelmatrix(q, k, Image)
% [H] = cupcoe_sphericalhankelmatrix(q, k, Image)
% -------------------------------------------------------------------------
% File:     couplingcoefficients.m
% Date:     21-09-18
% Author:   Lachlan Birnie
%
% Description:
%
%   Returns the spherical Hankel function in vector form for a given order
%   and set of image source radial distance properties.
%
%   - [H] = [1 by 8(2D+1)^3]
%
% Inputs:
%
%   q       Scalar integer, order of Hankel function terms.
%           
%   k       Scalar, wave number (frequency) in m^-1.
%           @TODO: allow function to accept vector k and output a matrix
%                  of spherical Hankel elements, [length(k) by 8(2D+1)^3]. 
%
%   Image   Structure with fields .r, the (image) radius aruguments.
%
% Outputs:
%
%   H       [1 by 8(2D+1)^3] matrix, where each element is:
%
%           H(1,b) = h_q (k* Image.r(b))
%
% Change Log:
%
%   26-09-18    - H is now a row vector, was a matrix.
%               - H no longer takes in Harm structure.
%
%% ------------------------------------------------------------------------

    % solve spherical harmonic equation for all images:
    H = sqrt(pi/2) .* 1./sqrt(k.*Image.r) .* besselh(q+.5, k.*Image.r);
       
end
% end    [H] = cupcoe_sphericalhankelmatrix(q, k, Image)

function [W] = cupcoe_wignermatrix(q, Harm, Image)
% [W] = cupcoe_wignermatrix(q, Harm, Image)
% -------------------------------------------------------------------------
% Function: wignermatrix
% File:     couplingcoefficients.m
% Author:   Lachlan Birnie
% Date:     27-09-18
%
% Desctiption:
%
%   Returns the 'Wigner terms' of the Hankel addition theorem for a given 
%   'q', where q = 0 -> N+V.
%
%   The 'Wigner terms' are: W = W1*W2*Z
%
%   where,  W1 = / n v q \
%                \ 0 0 0 / 
%
%           W2 = / n   v   q    \
%                \ m' -u (u-m') /
%
%           Z = sqrt( (2n+1) (2v+1) (2q+1) / (4pi) )
%
%  and,     m' = (-1)^(p1+p2) * m
%
% Inputs:
%
%   q       q^th index / integer in the Hankel addition theorem summation,
%           where q goes 0 -> N+V.
%
%   Harm    Struct with fields of .n .m .v .u for 'Wigner terms' arguments.
%
%   Image   Struct with fields .p1 .p2 and .num to get m'
%
% Ouputs:
%
%   W       [(V+1)^2 (N+1)^2 by 8(2D+1)^3] matrix of 'Wigner terms' for 
%           all harmonic + image pairs, for a given 'q'. 
%
%% ------------------------------------------------------------------------

    % solve Z for the harmonic indexes:
    Z = sqrt( ((2.*Harm.n+1) .* (2.*Harm.v+1) .* (2*q+1)) ./ (4*pi) );
    
    % find indexes of positive images, m' = m:
    i_pimg = (-1).^(Image.p1 + Image.p2) == 1;
    % find indexes of negative images, ,m' = -m:
    i_nimg = (-1).^(Image.p1 + Image.p2) == -1;
    
    % solve W1:
    W1 = cupcoe_wigner3jvector(Harm.n, Harm.v, ones(length(Harm.n),1).*q, ...
            zeros(Harm.num,1),zeros(Harm.num,1),zeros(Harm.num,1));
    
    % solve W2 when m' = m, (positive images):
    posW2 = cupcoe_wigner3jvector(Harm.n, Harm.v, ones(Harm.num,1).*q, ...
                           Harm.m, -Harm.u, (Harm.u-Harm.m) );
    
    % solve W2 when m' = -m, (negative images):
    minW2 = cupcoe_wigner3jvector(Harm.n, Harm.v, ones(Harm.num,1).*q, ...
                           -Harm.m, -Harm.u, (Harm.u--Harm.m) );
                       
    % @NOTE: used ones().*q to expand q to a column vector that matches
    %        the Harm. indexes, required for 'wigner3jvector' function. 
    
    % define final 'Wigner terms' matrix (invalid arguments remain zero):
    W2 = zeros(Harm.num, Image.num);
    
    % add positive and negative image W2's together:
    W2 = W2 + i_pimg .* posW2;
    W2 = W2 + i_nimg .* minW2;
    
    % solve 'wigner terms':
    W = W1 .* W2 .* Z;
    
end
% end    [W] = cupcoe_wignermatrix(q, Harm, Image)

function [W] = cupcoe_wigner3jvector(j1, j2, j3, m1, m2, m3)
% [W] = cupcoe_wigner3jvector(j1, j2, j3, m1, m2, m3)
% -------------------------------------------------------------------------
% Project   : wigner3jvector
% File      : couplingcoefficients.m
% Author    : Lachlan Birnie
% Creation  : 27-09-18
% Version   : v7.20 (30-09-19)
%
% Description:
%
%   Compute Wigner 3j symbol using Racah formula on the elements of input
%   vector terms.
%
% Inputs: 
%
%   / j1 j2 j3 \
%   \ m1 m2 m3 /
%
%       - where all terms are column vectors of the same length.
%
% Outputs:
%
%   W   Column vector of Wigner 3j symbols corresponding to the input
%       vectors.
%
%       W(i) = / j1(i) j2(i) j3(i) \
%              \ m1(i) m2(i) m3(i) /
%
% Acknowledgements:
%
%   Inspired by Wigner3j.m by Kobi Kraus, Technion, 25-06-08
%   - avaliable at: https://au.mathworks.com/matlabcentral/fileexchange
%                   /20619-wigner3j-symbol
%
% History
%
%   v7.20 30-09-19  - Fixed input validation criteria. The mistake had no
%                     effect on the coupling coefficient calculation, 
%                     but if this wigner function was used for another 
%                     project then the error may currupt the results.
%
% -------------------------------------------------------------------------

    % Find invalid inputs from the vectors:
    i_invalid = ( (j3 > (j1 + j2)) ...          % j3 out of interval
                | (j3 < abs(j1 - j2)) ...       % j3 out of interval
                | ((m1 + m2 + m3) ~= 0) ...     % non-conserving agular momentum 
                | (abs(m1) > j1) ...            % m is larger than j
                | (abs(m2) > j2) ...            % m is larger than j
                | (abs(m3) > j3) ...            % m is larger than j
                | ((~any([m1,m2,m3],2)) & (mod(j1+j2+j3,2)==1)) ... % m1 = m2 = m3 = 0 & j1 + j2 + j3 is odd
                );
    
    % find valid inputs:
    i_valid = ~i_invalid;
    n_valid = sum(i_valid);
    
    % return 0 if no valid inputs:
    if (n_valid == 0)
        W = zeros(size(j1));
        return;
    end
        
    % Select valid inputs:
    vj1 = j1(i_valid);
    vj2 = j2(i_valid);
    vj3 = j3(i_valid);
    vm1 = m1(i_valid);
    vm2 = m2(i_valid);
    vm3 = m3(i_valid);
    
    % compute terms:
    t1 = vj2 - vm1 - vj3;
    t2 = vj1 + vm2 - vj3;
    t3 = vj1 + vj2 - vj3;
    t4 = vj1 - vm1;
    t5 = vj2 + vm2;
    
    tmin = max( 0,  max( t1, t2 ) );
    tmax = min( t3, min( t4, t5 ) );
    
    % find largets summation:
    n_t =  max(tmax-tmin)+1;
    
    % fill in summation term matrix, pad with NaN.
    t = zeros(n_valid, n_t);
    for i = 1:n_valid
        t(i,:) = [(tmin(i):tmax(i)),nan(1,n_t-length(tmin(i):tmax(i)))];
    end
    
    % more terms?
    x = zeros(n_valid,n_t,6);
    x(:,:,1) = t;
    x(:,:,2) = t-t1;
    x(:,:,3) = t-t2;
    x(:,:,4) = t3-t;
    x(:,:,5) = t4-t;
    x(:,:,6) = t5-t;
    
    x2 = [ vj1+vj2+vj3+1 ...
           ,vj1+vj2-vj3 ...
           ,vj1-vj2+vj3 ...
           ,-vj1+vj2+vj3 ...
           ,vj1+vm1 ...
           ,vj1-vm1 ...
           ,vj2+vm2 ...
           ,vj2-vm2 ...
           ,vj3+vm3 ...
           ,vj3-vm3];
    
    % solve everything summation terms:
    sterm = (-1).^t .* exp( sum(-gammaln( x+1 ), 3) ...
                            + sum([-1,ones(1,9)].*gammaln( x2 +1 ), 2) ...
                            .* 0.5 );
                        
    % change NaNs to zero:
    sterm(isnan(sterm)) = 0;
    
    % sum over t:
    w_valid = sum(sterm, 2) .* (-1).^(vj1-vj2-vm3);
    
    % map back to full size:
    W = zeros(length(j1),1);
    W(i_valid) = w_valid;
       
end
% end    [W] = cupcoe_wigner3jvector(j1, j2, j3, m1, m2, m3)