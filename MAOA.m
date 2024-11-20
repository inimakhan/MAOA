

%__________________________________________________________________ %
%                    Archive-based Multi-Objective                  %
%               Arithmetic Optimization Algorithm (MAOA)           %
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



%% Function Details
function [Archive_costs]=MAOA(M_Iter,Archive_size,xNew_num,nVar,method,m)


if method==3

    TestProblem=sprintf('P%d',m);

    fobj = Ptest(TestProblem);

    xrange  = xboundaryP(TestProblem);
    nVar=max(size(xrange));
    % Lower bound and upper bound
    lb=xrange(:,1)';
    ub=xrange(:,2)';

end



%% Initialize the positions of solution

C_Iter=1;
Mu=0.1;
alpha=5;

MOP_Max=1;
MOP_Min=0.2;

Alpha=0.1;  % Grid Inflation Parameter
nGrid=30;   % Number of Grids per each Dimension
beta=4;     % Leader Selection Pressure Parameter
gamma=2;

Xnew=CreateEmptyParticle(xNew_num);


for i=1:xNew_num
    Xnew(i).Velocity=0;
    Xnew(i).Position=zeros(1,nVar);

    for j=1:nVar
        Xnew(i).Position(1,j)=unifrnd(lb(j),ub(j),1);
    end

    Xnew(i).Cost=fobj(Xnew(i).Position');
    Xnew(i).Best.Position=Xnew(i).Position;
    Xnew(i).Best.Cost=Xnew(i).Cost;
end

Xnew=DetermineDominations(Xnew);

Archive=GetNonDominatedParticles(Xnew);

Archive_costs=GetCosts(Archive);
G=CreateHypercubes(Archive_costs,nGrid,Alpha);

for i=1:numel(Archive)
    [Archive(i).GridIndex Archive(i).GridSubIndex]=GetGridIndex(Archive(i),G);
end

%% Main Loop

while C_Iter<M_Iter+1                        % Main loop
    MOP=1-((C_Iter)^(1/alpha)/(M_Iter)^(1/alpha));
    MOA=MOP_Min+C_Iter*((MOP_Max-MOP_Min)/M_Iter);
    %Upate the Position of solutions
    for i=1:xNew_num
        Leader=SelectLeader(Archive,beta);% if each of the ub and lb has a just value
        for j=1:nVar
            r1=rand();
            if (size(lb,2)==1)
                if r1<MOA
                    r2=rand();
                    if r2>0.5
                        Xnew(i).Position(1,j)=(Leader(1).Position(1,j)/(MOP+eps)*((ub(j)-lb(j))*Mu+lb(j)));
                    else
                        Xnew(i).Position(1,j)=(Leader(1).Position(1,j)*MOP*((ub(j)-lb(j))*Mu+lb(j)));
                    end
                else
                    r3=rand();
                    if r3>0.5
                        Xnew(i).Position(1,j)=(Leader(1).Position(1,j)-MOP*((ub(j)-lb(j))*Mu+lb(j)));
                    else
                        Xnew(i).Position(1,j)=(Leader(1).Position(1,j)+MOP*((ub(j)-lb(j))*Mu+lb(j)));
                    end
                end
            end


            if (size(lb,2)~=1)                          % if each of the ub and lb has more than one value
                r1=rand();
                if r1<MOA
                    r2=rand();
                    if r2>0.5
                        Xnew(i).Position(1,j)=(Leader(1).Position(1,j)/(MOP+eps)*((ub(j)-lb(j))*Mu+lb(j)));
                    else
                        Xnew(i).Position(1,j)=(Leader(1).Position(1,j)*MOP*((ub(j)-lb(j))*Mu+lb(j)));
                    end
                else
                    r3=rand();
                    if r3>0.5
                        Xnew(i).Position(1,j)=((Leader(1).Position(1,j)-MOP*((ub(j)-lb(j))*Mu+lb(j))));
                    else
                        Xnew(i).Position(1,j)=((Leader(1).Position(1,j)+MOP*((ub(j)-lb(j))*Mu+lb(j))));
                    end
                end
            end
            Xnew(i).Position=min(max(Xnew(i).Position,lb),ub);
            Xnew(i).Cost=fobj(Xnew(i).Position');
        end
    end
    Xnew=DetermineDominations(Xnew);
    non_dominated_Xnew=GetNonDominatedParticles(Xnew);
    Archive=[Archive
             non_dominated_Xnew];
    Archive=DetermineDominations(Archive);
    Archive=GetNonDominatedParticles(Archive);

    for i=1:numel(Archive)
        [Archive(i).GridIndex Archive(i).GridSubIndex]=GetGridIndex(Archive(i),G);
    end

    if numel(Archive)>Archive_size
        EXTRA=numel(Archive)-Archive_size;
        Archive=DeleteFromRep(Archive,EXTRA,gamma);
        Archive_costs=GetCosts(Archive);
        G=CreateHypercubes(Archive_costs,nGrid,Alpha);

    end
    % Display the iteration and best optimum obtained so far
    disp(['In iteration ' num2str(C_Iter) ': Number of solutions in the archive = ' num2str(numel(Archive))]);
    save results

    % Results
    Archive_costs=GetCosts(Archive);

    C_Iter=C_Iter+1;  % incremental iteration

end

end
