% Simulating 2-dimensional system with chemical reactions of type A+A -> 0
% In this script, you can choose number of particles (numP), the system
% length (L) and the number of experiments (numExperiments).

% By Safi Rafie-Zinedine

% The number of particles
numP = 20;
% The system length
L = 10;
% The number of dimensions
dim = 2;
% The volume of the system
volume = L^dim;
% The number of experiments
numExperiments = 1;

% Initializing the lattice 2D
Lattice_2D = zeros(L,L);
% Storing final lattices after each experiment
all_Lattice_2D = cell(1,numExperiments);

% The time step and time step for each particle
time_step = 1;
time_step_particle = time_step/numP;
% The start-time and the stop-time
t = 0;
time_stop = 300;
% Storing final time after each experiment
all_t = zeros(1,numExperiments);
% Storing the number of remaining particles after for each experiment
all_num = zeros(1,numExperiments);

% The density of particles
density = numP/volume;
% The density of particles decay with time
density_decay = zeros(1,time_stop+1);
% Storing the density after each experiment
all_density = cell(1,numExperiments);

% Choosing to occupy each site in lattice with probability num/volume
% for two different particle types: 
% num/volume = 2/randi_max  -> randi_max = 2*volume/num
randi_max = round(2*volume/numP);

% Looping through all lattice sites and fill them randomly with equal
% probability
total_numP = 0;
while total_numP~=numP
    for i=1:L
        for j=1:L
            if ((Lattice_2D(i,j) ~= 1) || (Lattice_2D(i,j) ~= 2))
               Lattice_2D(i,j) = randi(randi_max);
               [iA,~] = ind2sub(size(Lattice_2D),find( Lattice_2D==1));
               [iB,~] = ind2sub(size(Lattice_2D),find( Lattice_2D==2));
               % The number of A particles
               numP_A = length(iA);
               % The number of B particles
               numP_B = length(iB);
               % The total number
               total_numP = numP_A + numP_B;
            end
            if total_numP == 20, break; end
        end
        if total_numP == 20, break; end
    end
end

% All ones particles (A)
isOne = Lattice_2D==1;
[iA,jA] = ind2sub(size(Lattice_2D),find( Lattice_2D==1)); 
% All twos particle (B)
isTwo = Lattice_2D==2;
[iB,jB] = ind2sub(size(Lattice_2D),find( Lattice_2D==2));
% All particles in the lattice
Lattice_2D = isOne + 2*isTwo;

% Creating the particle in the lattice
% First subplot to see the initial state of the particle in the lattice
create_particles(1,L,iA,jA,iB,jB);
% Second subplot to update the particle state in the lattice
C = create_particles(2,L,iA,jA,iB,jB);
% Third subplot to plot the particle density with time.
D = create_density(time_stop,t,density);

% The number of experiments
for E = 1:numExperiments
    % Re-initializing variables for each experiment
    updated_Lattice_2D = Lattice_2D;
    updated_numP = numP;
    t = 0;    
    for T = 0:time_step:time_stop
    % Choosing random particle in the lattice
    particle = [randi(L), randi(L)];
    while updated_Lattice_2D(particle(1,1), particle(1,2))==0
        particle(1,1) = randi(L);
        particle(1,2) = randi(L);
    end

    % Determining if the chosen particle of type A or type B
    chosen_particle = updated_Lattice_2D(particle(1,1), particle(1,2));
    if chosen_particle == 1
        anti_chosen_particle = 2;
    else
        anti_chosen_particle = 1;
    end
    
    % Choosing random nearest neighbor site
    neighbor = randi(4);
    switch neighbor
        case 1
            site_near = particle + [0, 1];  % particle_up
        case 2
            site_near = particle + [0, -1]; % particle_down
        case 3
            site_near = particle + [1, 0]; % particle_right
        case 4
            site_near = particle + [-1, 0]; % particle_left
    end

    % Ensuring periodic boundaries of the system
    if site_near(1,1) <= 0
        site_near(1,1) = L;
    elseif site_near(1,1) > L
        site_near(1,1) = 1;
    elseif site_near(1,2) > L
        site_near(1,2) = 1;
    elseif site_near(1,2) <= 0
        site_near(1,2) = L;
    end

    % Killing particles if they meet in the same site, otherwise updating
    % their sites
    if updated_Lattice_2D(site_near(1,1), site_near(1,2)) == chosen_particle     
       updated_Lattice_2D(site_near(1,1), site_near(1,2)) = 0;
       updated_Lattice_2D(particle(1,1), particle(1,2)) = 0;
       updated_numP = updated_numP - 2;
    elseif updated_Lattice_2D(site_near(1,1), site_near(1,2)) == anti_chosen_particle
       updated_Lattice_2D(site_near(1,1), site_near(1,2)) = 0;
       updated_Lattice_2D(particle(1,1), particle(1,2)) = 0;
       updated_numP = updated_numP - 2;        
    else
       updated_Lattice_2D(site_near(1,1), site_near(1,2)) =...
                updated_Lattice_2D(particle(1,1), particle(1,2));
       updated_Lattice_2D(particle(1,1), particle(1,2)) = 0;
    end
    
    % Storing and ploting particle density
    t = t + time_step;
    density_decay(1,T+1) = updated_numP/volume;
    update_density(time_stop,density_decay,D);
    
    % Storing and ploting particles
    [iExp1A,jExp1A] = ind2sub(size(updated_Lattice_2D),find(updated_Lattice_2D==1));
    [iExp1B,jExp1B] = ind2sub(size(updated_Lattice_2D),find(updated_Lattice_2D==2));
    update_particles(iExp1A,jExp1A,iExp1B,jExp1B,C);
    
    % Breaking the loop if all particle have been killed
    if updated_numP==0, break; end
    
    end
    
    % Storing and ploting particles and particle density after each experiment.
    all_num(1,E) = updated_numP;
    all_t(1,E) = t;
    all_Lattice_2D{E} = updated_Lattice_2D;
    all_density{E} = density_decay;    
end

% Creating and ploting particles in lattice
function C = create_particles(num_fig,L,iExp1A,jExp1A,iExp1B,jExp1B)
   % Create the appropriate figure with axes with grid
   subplot(2,2,num_fig)
   set(gcf,'Position', [300,350,1400,950]);
   title('(A)');
   ax = gca;
   ax.YLim = [0 L+1];
   xlabel('X  (\mum)');
   ax.XLim = [0 L+1];
   ylabel('Y  (\mum)');
   grid on
   
   % Ploting
   C.particle_A = line(iExp1A,jExp1A,'LineStyle','none','marker','.',...
       'Color','b','markersize',40);
   C.particle_B = line(iExp1B,jExp1B,'LineStyle','none','marker','.',...
       'Color','b','markersize',40);
end

% Updating and reploting the particles
function update_particles(iExp1A,jExp1A,iExp1B,jExp1B,C)         
   set(C.particle_A,'XData',iExp1A,'YData',jExp1A);
   set(C.particle_B,'XData',iExp1B,'YData',jExp1B);
   drawnow       
end

% Ploting particle density
function D = create_density(time_stop,t,density)
   % Create the appropriate figure with axes with grid
   subplot(2,2,[3 4])
   title('(C)');
   ax = gca;
   ax.YLim = [0 density];
   xlabel('Time (sec)');
   ax.XLim = [0 time_stop];
   ylabel('Particle Density  (number/volume)');
   hold(ax,'on');
   grid on

   % Ploting
   D.densityExp = plot(t,density,'b.','MarkerSize',8);
end

% Reploting particle density
function update_density(time_stop,density_decay,D)         
   set(D.densityExp,'XData',0:time_stop,'YData',density_decay);
   drawnow       
end