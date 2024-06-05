# GeneralizedISMCouplingCoefficients
Spherical harmonic region-to-region image source method coupling coefficients


    I M A G E - S O U R C E - C O U L P I N G - C O E F F I C I E N T S 
                                ( +CupCoe )


Package:    +CupCoe
Creation:   20-09-18
Version:    7.30 (29-11-19)
Authors:    Lachlan Birnie
            Prasanga N. Samarasinghe 
            Thushara D. Abhayapala
            Yan Lu
            Glenn Dickins
            Hanchi Chen

DESCRIPTION

    MATLAB package to calculate the image source method based coupling 
    coefficients (alpha_{\nu\mu}^{nm}(k)), of the spherical harmonic 
    based region-to-region room transfer function.

    As presented in Eq. (13) of:

        Samarasinghe, P.N., Abhayapala, T.D., Lu, Y., Chen, H. 
        and Dickins, G., 2018. Spherical harmonics based generalized 
        image source method for simulating room acoustics. 
        The Journal of the Acoustical Society of America, 
        144(3), pp.1381-1391.

___________________________________________________________________________
GETTING - STARTED

    All required MATLAB files are contained within the "+CupCoe" folder.
    - Copy the +CupCoe folder into the same file (working directory) as
      your .m script/function. 

    
    To generate a set of (alpha) coupling coefficients:
    - For a given frequency defined by 'k', 
    - a given room (size, reflections, image depth) defined by 'Room', 
    - a given source region (origin, radius, order) defined by 'Source',
    - a given receiver (origin, radius, order) region defined by 'Region',
    - Call:

        [A] = CupCoe.couplingcoefficients(k, Room, Source, Receiver);

        - "CupCoe." informs MATLAB that the couplingcoefficients()
          function is contained in the +CupCoe folder. 
        - A is a [(V+1)^2 by (N+1)^2] matrix of coupling coefficeints.


    Alternatively, you can call:

        [A] = CupCoe.load_couplingcoefficients(k, Room, Source, Receiver);

        - This wrapper function will save the [A] matrix to a .mat file
          in \+CupCoe\database\ when it is first calculated.
        - Calling the load_couplingcoefficients() function a second time
          with the same inputs will load [A] from the database, instead 
          of recalculating them again.

___________________________________________________________________________
U S I N G - C O U P L I N G - C O E F F I C I E N T S

    The matrix of alpha coupling coefficients is denoted as [A]:
    - [A] is a 2D matrix of complex elements and size [(V+1)^2 by (N+1)^2]
    - Each element of [A] is a coupling coefficient:


                    | a_(00)^(00)    ...    a_(00)^(NN) |
                    |      .      .               .     |
        [A](k) =    |      .    a_(v'u')^(n'm')   .     |
                    |      .              .       .     |
                    | a_(VV)^(00)    ...    a_(VV)^(NN) |

    
    - A column of [A] is denoted as (for a unit source mode 00):

        A(:,1) = [a_(00)^(00) a_(1-1)^(00) a(10)^(00) a(1+1)^(00) ... ]^T

    - A row of [A] is denoted as (for a unit receiver mode 00):
        
        A(1,:) = [a_(00)^(00) a_(00)^(1-1) a_(00)^(10) a_(00)^(1+1) ... ]

    - To get a specific coupling coefficient:

        a_(vu)^(nm) = A(v^2+v+u+1, n^2+n+m+1);

___________________________________________________________________________
FILES

    +CupCoe\couplingcoefficeints.m
        - This is the main function which calculates the (alpha) coupling 
          coefficients for a given system. 

    +CupCoe\load_couplingcoefficients.m
        - This is a wrapper function for couplingcoefficients.m, which 
          handles the saving and loading of [A] matrices. 
        - When coupling coefficeints are calculated for the first time,
          they are saved to the +cupcoe\database folder with a unique
          .mat file name. 
        - Subsequent calls with the same inputs will load the coupling 
          coefficients from the database folder instead of re-calculating.

    +CupCoe\database\
        - Folder containing [A] matrices as .mat files. 
        - Folders creation is handled by load_ wrapper function. 
        - Example .mat file name: 
            cupcoe_k18p5_L463D3_B0p30-0p40-0p50-0p60-0p70
            -0p80_OS1p01p51p5RS1p0N19_OR2p53p01p5RR0p0V3.mat

    +CupCoe\Smarasinghe_etal_Spherical-harmonic-based-generalized-image-
    source-method_JASA_2018.pdf
        - Copy of paper detailing the coupling coefficient equation. 

    +CupCoe\readme.txt
        - This file.

___________________________________________________________________________
DEPENDENCIES

    MATLAB  :   The +CupCoe Package has been tested with
                MATLAB R2018a
                MATLAB R2017a

                (MATLAB 2016a is known to not work)


    The couplingcoefficient() function depends on 4 sub-functions 
    contained within the +CupCoe\couplingcoefficients.m file:

        cupcoe_sphericalharmonicmatrix(q, Harm, Image)
        cupcoe_sphericalhankelmatrix(q, k, Image)
        cupcoe_wignermatrix(q, Harm, Image)
        cupcoe_wigner3jvector(j1, j2, j3, m1, m2, m3)

___________________________________________________________________________
E X A M P L E 

    % copy paste this code - - - - - - - - - - - - - - - - - - - - - - - -
        k = 18.3;                           % wave number (frequency)

        Room.size = [4 6 3];                % (Lx, Ly, Lz) in meters
        Room.b = [.9 .7 .8 .95 .78 .81];    % reflection coefficients
        Room.D = 3;                         % image depth

        Source.Os = [1 1 1];                % (x,y,z) wrt bot-front-left 
        Source.Rs = 1;                      % radial size in meters
        Source.N = ceil(k*Source.Rs);       % truncation order

        Receiver.Or = [2.1 2.1 2.1];        % (x,y,z) wrt bot-front-left 
        Receiver.Rr = 0.042;                % radial size in meters
        Receiver.V = 3;                     % truncation order

        [A] = CupCoe.couplingcoefficients(k, Room, Source, Receiver);
    %end - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

___________________________________________________________________________
HISTORY

v2.00 : 31-10-18 : Lachlan Birnie : Release MATLAB package version (v2.00)
                                  : couplingcoefficients.m          v7.01
                                  : load_couplingcoefficients.m     v1.00
        02-11-18 : Lachlan Birnie : Renamed package to +CupCoe

v2.01 : 05-11-18 : Lachlan Birnie : Fixed some +CupCoe notation typos.
                                  : load_couplingcoefficients.m     v1.01

v2.02 : 06-11-18 : Lachlan Birnie : Documentation updates.
                                  : couplingcoefficients.m          v7.02

v2.10 : 30-09-19 : Lachlan Birnie : Fixed Wigner function input validation.
                                  : couplingcoefficeints.m          v7.20

-- match version number with couplingcoefficients.m -- 

v7.30 : 29-11-19 : Lachlan Birnie : Fixed bug when sources on z-axis. 

