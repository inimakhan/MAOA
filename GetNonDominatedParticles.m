

%__________________________________________________________________ %
%                    Archive-based Multi-Objective                  %
%               Arithmetic Optimization Algorithm (MAOA)            %
%                                                                   %
%                                                                   %
%                  Developed in MATLAB R2022a (MacOs)               %
%                                                                   %
%                     Author and programmer                         %
%                ---------------------------------                  %
%                      Nima Khodadadi (ʘ‿ʘ)                         %
%                             e-Mail                                %
%                ---------------------------------                  %
%                         nkhod002@fiu.edu                          %
%                                                                   %
%                            Homepage                               %
%                ---------------------------------                  %
%                    https://nimakhodadadi.com                      %
%                                                                   %
%                                                                   %
%                                                                   %
%                                                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% ----------------------------------------------------------------------- %









function nd_pop=GetNonDominatedParticles(pop)

    ND=~[pop.Dominated];
    
    nd_pop=pop(ND);

end