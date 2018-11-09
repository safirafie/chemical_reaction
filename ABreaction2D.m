% Simulating 2-dimensional system with chemical reactions of type A+B -> 0
% In this script, you can choose number of particles (numP), the system
% length (L) and the number of experiments (numExperiments).

% By Safi Rafie-Zinedine

% Number of particles
numP = 20;
% The system length
L = 10;
% Number of dimensions
dim = 2;
% Volume of the system
volume = L^dim;
% Number of experiments
numExperiments = 1;

% Initialize the lattice 2D
Lattice_2D = zeros(L,L);
% Final lattices after numExperiments
all_Lattice_2D = cell(1,numExperiments);

% The time step
time_step = 1;
time_step_particle = time_step/numP;
% The start-time and the stop-time
t = 0;
time_stop = 300;
all_t = zeros(1,numExperiments);
% Collecting remaining particle number for each experiment
all_num = zeros(1,numExperiments);

% density of particles
density = numP/volume;
% The density of particles decay with time
density_decay = zeros(1,time_stop+1);
% Resulting density after each experiment
all_density = cell(1,numExperiments);

numP_A = 0;
% The density of A particles decay with time
density_decay_A = zeros(1,time_stop+1);
% Resulting A particle density after each experiment
all_density_A = cell(1,numExperiments);

numP_B = 0;
% The density of B particles decay with time
density_decay_B = zeros(1,time_stop+1);
% Resulting B particle density after each experiment
all_density_B = cell(1,numExperiments);

% Choose to occupy that site with probability num/volume
% for two different particle types.
% num/volume = 2/randi_max  -> randi_max = 2*volume/num
randi_max = round(2*volume/numP);

total_numP = 0;
while total_numP~=numP
    for i=1:L
        for j=1:L
            if ((Lattice_2D(i,j) ~= 1) || (Lattice_2D(i,j) ~= 2))
               Lattice_2D(i,j) = randi(randi_max);
               [iA,~] = ind2sub(size(Lattice_2D),find( Lattice_2D==1));
               [iB,~] = ind2sub(size(Lattice_2D),find( Lattice_2D==2));
               % Number of A particles
               numP_A = length(iA);
               % Number of B particles
               numP_B = length(iB);
               % Total number
               total_numP = numP_A + numP_B;
            end
            if total_numP == 20, break; end
        end
        if total_numP == 20, break; end
    end
end

% all ones particles (A)
isOne = Lattice_2D==1;
[iA,jA] = ind2sub(size(Lattice_2D),find( Lattice_2D==1)); 
% all twos particle (B)
isTwo = Lattice_2D==2;
[iB,jB] = ind2sub(size(Lattice_2D),find( Lattice_2D==2));
% all particles in the lattice
Lattice_2D = isOne + 2*isTwo;
% density of A particles and B particles
density_A = numP_A/volume;
density_B = numP_B/volume;

create_particles(1,L,iA,jA,iB, jB);
C = create_particles(2,L,iA,jA,iB, jB);
D = create_density(time_stop,t,density,density_A,density_B);

% The number of experiments
for E = 1:numExperiments
    % Initialize variable for each experiment
    updated_Lattice_2D = Lattice_2D;
    updated_numP = numP;
    updated_numP_A = numP_A;
    updated_numP_B = numP_B;
    t = 0;
    for T = 0:time_step:time_stop
    % Choosing random particle
    particle = [randi(L), randi(L)];
    while updated_Lattice_2D(particle(1,1), particle(1,2))==0
        particle(1,1) = randi(L);
        particle(1,2) = randi(L);
    end

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

    % Ensure periodic boundary
    if site_near(1,1) <= 0
        site_near(1,1) = L;
    elseif site_near(1,1) > L
        site_near(1,1) = 1;
    elseif site_near(1,2) > L
        site_near(1,2) = 1;
    elseif site_near(1,2) <= 0
        site_near(1,2) = L;
    end

    % Kill particles if double occupied site
    if updated_Lattice_2D(site_near(1,1), site_near(1,2)) == chosen_particle
       updated_Lattice_2D(site_near(1,1), site_near(1,2)) = chosen_particle;
       updated_Lattice_2D(particle(1,1), particle(1,2)) = chosen_particle;      
    elseif updated_Lattice_2D(site_near(1,1), site_near(1,2)) == anti_chosen_particle
       updated_Lattice_2D(site_near(1,1), site_near(1,2)) = 0;
       updated_Lattice_2D(particle(1,1), particle(1,2)) = 0;
       updated_numP = updated_numP - 2;
       updated_numP_A = updated_numP_A - 1;
       updated_numP_B = updated_numP_B - 1;
    else
       updated_Lattice_2D(site_near(1,1), site_near(1,2)) =...
                updated_Lattice_2D(particle(1,1), particle(1,2));
       updated_Lattice_2D(particle(1,1), particle(1,2)) = 0;
    end
    
    t = t + time_step;
    density_decay(1,T+1) = updated_numP/volume;
    density_decay_A(1,T+1) = updated_numP_A/volume;
    density_decay_B(1,T+1) = updated_numP_B/volume;
    update_density(time_stop,density_decay_A,density_decay_B,D);
    
    [iExp1A,jExp1A] = ind2sub(size(updated_Lattice_2D),find(updated_Lattice_2D==1));
    [iExp1B,jExp1B] = ind2sub(size(updated_Lattice_2D),find(updated_Lattice_2D==2));
    update_particles(iExp1A,jExp1A,iExp1B,jExp1B,C);
    
    if updated_numP==0, break; end
    
    end
    all_num(1,E) = updated_numP;
    all_t(1,E) = t;
    all_Lattice_2D{E} = updated_Lattice_2D;
    all_density{E} = density_decay;
    all_density_A{E} = density_decay_A;
    all_density_B{E} = density_decay_B;
end

% Adding legend at the end.
subplot(2,2,[3 4]);
legend('Particle Density','Anit-particle Density');

% create and plot function
function C = create_particles(num_fig,L,iExp1A,jExp1A,iExp1B,jExp1B)
   % Create the appropriate figure with axes with grid
   subplot(2,2,num_fig)
   set(gcf,'Position', [300,350,1400,950]);
   title('(A)');
   ax = gca;
   ax.YLim = [0 L+1];
   ax.XLim = [0 L+1];
   xlabel('X  (\mum)');
   ylabel('Y  (\mum)');
   grid on
       
   C.particle_A = line(iExp1A,jExp1A,'LineStyle','none','marker','.',...
       'Color','b','markersize',40);
   C.particle_B = line(iExp1B,jExp1B,'LineStyle','none','marker','.',...
       'Color','r','markersize',40);
end

% Update and replot function(position,C)
function update_particles(iExp1A,jExp1A,iExp1B,jExp1B,C)         
   set(C.particle_A,'XData',iExp1A,'YData',jExp1A)
   set(C.particle_B,'XData',iExp1B,'YData',jExp1B)
   drawnow       
end

% create and plot density
function D = create_density(time_stop,t,density,density_A,density_B)
   % Create the appropriate figure with axes with grid
   subplot(2,2,[3 4]);
   title('(C)');
   ax = gca;
   ax.YLim = [0 density];
   ax.XLim = [0 time_stop];
   xlabel('Time (sec)');
   ylabel('Particle Density  (number/volume)');
   hold(ax,'on');
   grid on
       
   D.densityExp_A = plot(t,density_A,'b*','MarkerSize',8);
   D.densityExp_B = plot(t,density_B,'r.','MarkerSize',10);
end

% Update and replot density
function update_density(time_stop,density_decay_A,density_decay_B,D)         
   set(D.densityExp_A,'XData',0:time_stop,'YData',density_decay_A);
   set(D.densityExp_B,'XData',0:time_stop,'YData',density_decay_B);
   drawnow       
end