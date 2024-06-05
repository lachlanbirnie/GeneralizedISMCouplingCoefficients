function [A] = load_couplingcoefficients(k, Room, Source, Receiver)
% [A] = load_couplingcoefficients(k, Room, Source, Receiver)
% =========================================================================
% File:         +cupcoe\load_couplingcoefficients.m
% Creation:     31-10-18
% Author:       Lachlan Birnie
% Version:      v1.01 (05-11-18)
%
% Description:
%
%   This is a wrapper function for couplingcoefficients.m. 
%  
%   load_couplingcoefficients() will attempt to load a set of coupling 
%   coefficients for a system, if they have already been calculated and
%   saved to a folder. If not, this function will calculate and save the 
%   coupling coefficients for the system.
%
% Inputs / Ouputs:
%
%   - See +cupcoe\couplingcoefficients.m
%
% Dependencies:
%
%   +cupcoe\couplingcoefficients.m
%
% Notes:
%
%   - The coupling coefficient .mat file names are only unique to 1 
%     decimal place, so minor changes to system parameters may load 
%     incorrect coupling coefficients.
%
% History
%
%   v1.01   05-11-18    - Fixed '+cupcoe' to '+CupCoe' in CUPCOE_FOLDER_DIR
%                       - Fixed 'A = CupCoe.couplingcoefficients()' 
%
% =========================================================================

    % Coupling coefficient database name (can edit):        
    CUPCOE_FOLDER_DIR = [pwd '\+CupCoe\'];
    CUPCOE_FOLDER_NAME = 'database';
    
    cupcoeDatabase = [CUPCOE_FOLDER_DIR '\' CUPCOE_FOLDER_NAME];
    
    % Coupling coefficient .mat name:
    cupcoeName = ...
    [strrep( ...
    [sprintf('cupcoe_k%.1f',k)...
    ,sprintf('_L%.1f%.1f%.1fD%i',Room.size,Room.D) ...
    ,sprintf('_B%.2f-%.2f-%.2f-%.2f-%.2f-%.2f', ...
        Room.b(1),Room.b(2),Room.b(3),Room.b(4),Room.b(5),Room.b(6))...
    ,sprintf('_OS%.1f%.1f%.1fRS%.1fN%i',Source.Os,Source.Rs,Source.N)...
    ,sprintf('_OR%.1f%.1f%.1fRR%.1fV%i', ...
        Receiver.Or,Receiver.Rr,Receiver.V)...
    ],'.','p') ...
    ,'.mat'...
    ];

    % Maintain coupling coefficient database:
    if (exist(cupcoeDatabase, 'dir') ~= 7)
        % make new folder:
        fprintf('load_couplingcoefficients(): making new folder');
        mkdir(cupcoeDatabase);
    end
    
    % Check if coupling coefficients are in database:
    if (exist([cupcoeDatabase '\' cupcoeName], 'file') == 2)
        fprintf('load_couplingcoefficients(): loading: %s\n', cupcoeName);
        temp = load([cupcoeDatabase '\' cupcoeName], 'A');
        A = temp.A;
    else
    % does not exist - calculate coefficients and save:
        A = CupCoe.couplingcoefficients(k, Room, Source, Receiver);        
        fprintf('load_couplingcoefficients(): saving: %s\n', cupcoeName);
        save([cupcoeDatabase '\' cupcoeName], 'A');
    end
            
end % [A] = load_couplingcoefficients(k, Room, Source, Receiver)