%% Generates gcode from boundaries
% pseudocode:
% for each boundary
%   go to start position with G0
%   currPos = start;
%   for segment 2 -> boundary(end)
%       G1 to segment E(calculated with dist function)
%       currPos = segment;
%       add to grid and deal with intersections

%% 
% boundaries -> cell array from a bwboundaries call. scaled to mm
% scale -> ratio from pixels to mm
% filename -> output
function gen_gcode(boundaries, filename)
    FILAMENT_DIA = 1.74; %mm
    EXTRUDE_MULT = 1; %initial extra material
    LAYER_HEIGHT = 0.15; %mm
    PART_HEIGHT = 15; % in mm, gets raised to factor of layer height
    NOZZLE_DIA = 0.35; %mm
    TRAVEL_SPEED = 3000;%6000;
    INIT_WORK_SPEED = 1200;%1800;
    WORK_SPEED = 1800;%3000;
    RETRACT_LEN = 3;%3; %mm
    RETRACT_SPEED = 6000; %mm/min
    INIT_Z_HEIGHT = 2.96; %acrylic height, change in start gcode as well
    
    system(sprintf('echo > %s', filename));
    system(sprintf('cat start.gcode > %s', filename));
    %system(sprintf('cat ellipse.gcode >> %s', filename));
    f = fopen(filename, 'a');
    fprintf(f, ';End start gcode and Ellipse base\n\n');
    currPos = [0, 0, 0];
    currExtrude = 0;%6574.61658;
    currZ = INIT_Z_HEIGHT;
 
    
    if RETRACT_LEN > 0
        fprintf(f, 'G1 F%d E%0.5f\n', RETRACT_SPEED, currExtrude - RETRACT_LEN);
        currExtrude = currExtrude - RETRACT_LEN;
    end  
        
    % Start layer loop here
    for z = 1:ceil(PART_HEIGHT/LAYER_HEIGHT)
        currZ = currZ + LAYER_HEIGHT;
        
        % For each boundary
        for i = 1:length(boundaries)
            boundary = boundaries{i};
            
            currPos = [boundary(1,1), boundary(1,2), currZ];
            fprintf(f,';Layer %d, boundary %d...\n', z, i);
            fprintf(f, 'G0 F%d X%0.3f Y%0.3f Z%0.3f\n', ...
                TRAVEL_SPEED, currPos(1), currPos(2), currPos(3));
            
            if RETRACT_LEN > 0
                fprintf(f, 'G1 F%d E%0.5f\n', RETRACT_SPEED, currExtrude + RETRACT_LEN);
                currExtrude = currExtrude + RETRACT_LEN;
            end
            for j = 2:length(boundary)
                nextPos = [boundary(j,1), boundary(j,2), currZ];
                euc_dist = dist([currPos(1:2)', nextPos(1:2)']);
                euc_dist = euc_dist(2,1);
                
                ext_len = euc_dist * EXTRUDE_MULT * (NOZZLE_DIA*LAYER_HEIGHT/(pi*(FILAMENT_DIA/2)^2));
                nextExtrude = currExtrude + ext_len;
                
                % set feedrate at the start of printing boundary
                if j == 2
                    % change speed for first layer
                    if z == 1
                        fprintf(f, 'G1 F%d X%0.3f Y%0.3f E%0.5f\n', ...
                            INIT_WORK_SPEED, nextPos(1), nextPos(2), nextExtrude);
                    else 
                        fprintf(f, 'G1 F%d X%0.3f Y%0.3f E%0.5f\n', ...
                            WORK_SPEED, nextPos(1), nextPos(2), nextExtrude);
                    end
                else
                    fprintf(f, 'G1 X%0.3f Y%0.3f E%0.5f\n', ...
                        nextPos(1), nextPos(2), nextExtrude);
                end
                
                currPos = nextPos;
                currExtrude = nextExtrude;
            end
            % Retract at the end of a boundary
            if RETRACT_LEN > 0
                fprintf(f, 'G1 F%d E%0.5f\n', RETRACT_SPEED, currExtrude - RETRACT_LEN);
                currExtrude = currExtrude - RETRACT_LEN;
            end
        end
        fprintf(1,'%d...', z);
        EXTRUDE_MULT = 1;
    end
    fprintf(1,'\n');
    fprintf(f, ';End print \n\n');
    system(sprintf('cat end.gcode >> %s', filename));
    fclose(f);
end