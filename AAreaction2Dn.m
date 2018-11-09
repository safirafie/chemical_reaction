% Simulating 2-dimensional system with chemical reactions of type A+A -> 0
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
% Remaining particle after each experiment
all_numP = zeros(1,numExperiments);

% The time step
time_step = 1;
time_step_particle = time_step/numP;
% The start-time and the stop-time
t = 0;
time_stop = 300;
% Spending time in each experiment
all_t = zeros(1,numExperiments);

% density of particles
density = numP/volume;
% The density of particles decay with time
density_decay = zeros(1,time_stop+1);
% Resulting density after each experiment
all_density = cell(1,numExperiments);

% Resulting position after each experiment
all_position = cell(1,numExperiments);
% position = numP by 2 array of dim=2D position of N particles.
position = randi(L,numP,dim);
           
% Avoid double occupation
no_double_occupation = true;
count_double_occ = 0;
while no_double_occupation
     for i=1:numP-1
        for j=i+1:numP
           if (position(i,1)==position(j,1))...
                           && (position(i,2)==position(j,2))               
                 position(j,:) = randi(L,1,2);
                 position(i,:) = randi(L,1,2);
                 count_double_occ = count_double_occ + 1;
           end
        end
     end
     if count_double_occ==0
        disp('No double occupation')
        no_double_occupation = false;
     elseif count_double_occ~=0
        count_double_occ = 0;
     end
end

% C=create figure and particles then plot them.
iA = position(:,1);
jA = position(:,2);
create_particles(1,L,iA,jA);
C = create_particles(2,L,iA,jA);
D = create_density(time_stop,t,density);

% Number of experiments
for E = 1:numExperiments
    updated_position = position;
    updated_numP = numP;
    t=0;
     % T iterations
     for T = 0:time_step:time_stop
         % random number of velocities.
         velocity = round(randn(updated_numP,dim)/2);
         
         % Calculate the position. x= x0 + vt
         updated_position = updated_position + time_step*velocity;
         
         % Ensure periodic boundary  
           updated_position(updated_position > L) = 1;
           updated_position(updated_position <= 0) = L;         
         % Kill meeted particles
           for i=1:updated_numP-1
               for j=i+1:updated_numP
                   if (updated_position(i,1)==updated_position(j,1))...
                           && (updated_position(i,2)==updated_position(j,2))                     
                      updated_position(j,:)=[];
                      updated_position(i,:)=[];
                      updated_numP = updated_numP - 2;
                      break;
                   end
                   break;
               end
               if updated_numP==0, break; end
           end
           
           t = t + time_step;
           density_decay(1,T+1) = updated_numP/volume;
           update_density(time_stop,density_decay,D);
           
           % Update function for position
           iExp1A = updated_position(:,1);
           jExp1A = updated_position(:,2);
           update_particles(iExp1A,jExp1A,C);
           
           if updated_numP==0, break; end
     end
   all_numP(1,E) = updated_numP;
   all_t(1,E) = t;
   all_position{E} = updated_position;
   all_density{E} = density_decay;
end

% create and plot function
function C = create_particles(num_fig,L,iExp1A,jExp1A)
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
end

% Update and replot function(position,C)
function update_particles(iExp1A,jExp1A,C)         
   set(C.particle_A,'XData',iExp1A,'YData',jExp1A)
   drawnow       
end

% create and plot density
function D = create_density(time_stop,t,density)
   % Create the appropriate figure with axes with grid
   subplot(2,2,[3 4])
   title('(C)');
   ax = gca;
   ax.YLim = [0 density];
   ax.XLim = [0 time_stop];
   xlabel('Time (sec)');
   ylabel('Particle Density  (number/volume)');
   hold(ax,'on');
   grid on
       
   D.densityExp = plot(t,density,'b.','MarkerSize',8);
end

% Update and replot density
function update_density(time_stop,density_decay,D)         
   set(D.densityExp,'XData',0:time_stop,'YData',density_decay);
   drawnow       
end