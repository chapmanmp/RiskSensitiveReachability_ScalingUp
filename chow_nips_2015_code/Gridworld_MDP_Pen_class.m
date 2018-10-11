classdef Gridworld_MDP_Pen_class < Finite_MDP_class
    properties
        Nrow = [];
        Ncol = [];
        img = [];
        obstacles = [];
        nonobstacles = [];
        targetRow = [];
        targetCol = [];
        die_state = [];
        eps = 0;
    end
    methods(Static) % Static method used to define the class instance
        function [newrow, newcol,collision_index]= north(row,col,Nrow,Ncol,im)
            collision_index = 0;
            newrow = max(row-1,1);
            newcol = col;
            
            %%% if we hit an obstacle   collision_index = 1;
            if im(newrow, newcol) ==0
                collision_index = 1;
                newrow = Inf;
                newcol = Inf;
            end
        end
        
        function [newrow, newcol,collision_index]= south(row,col,Nrow,Ncol,im)
            collision_index = 0;
            newrow = min(row+1,Nrow);
            newcol = col;
            
            %%% if we hit an obstacle   collision_index = 1;
            if im(newrow, newcol) ==0
                collision_index = 1;
                newrow = Inf;
                newcol = Inf;
            end
        end
        
        function [newrow, newcol,collision_index]= east(row,col,Nrow,Ncol,im)
            collision_index = 0;
            newrow = row;
            newcol = min(col+1,Ncol);
            
            %%% if we hit an obstacle   collision_index = 1;
            if im(newrow, newcol) ==0
                collision_index = 1;
                newrow = Inf;
                newcol = Inf;
            end
        end
        
        function [newrow, newcol,collision_index]= west(row,col,Nrow,Ncol,im)
            collision_index = 0;
            newrow = row;
            newcol = max(col-1,1);
            
            %%% if we hit an obstacle   collision_index = 1;
            if im(newrow, newcol) ==0
                collision_index = 1;
                newrow = Inf;
                newcol = Inf;
            end
        end
        
        function [newrow,newcol,newcollision_index] = neigborhood_fn(row,col,Nrow,Ncol,im)
            % Generate neigbors of state (row,col)
            newrow = zeros([1 4]);
            newcol = zeros([1 4]);
            newcollision_index = zeros([1 4]);
            
            %north
            clear row_buf col_buf collision_index_buf
            [row_buf, col_buf, collision_index_buf]= Gridworld_MDP_Pen_class.north(row,col,Nrow,Ncol,im);
            newrow(1) = row_buf;
            newcol(1) = col_buf;
            newcollision_index(1) = collision_index_buf;
            
            %south
            clear row_buf col_buf collision_index_buf
            [row_buf, col_buf, collision_index_buf]= Gridworld_MDP_Pen_class.south(row,col,Nrow,Ncol,im);
            newrow(2) = row_buf;
            newcol(2) = col_buf;
            newcollision_index(2) = collision_index_buf;
            
            %east
            clear row_buf col_buf collision_index_buf
            [row_buf, col_buf, collision_index_buf]= Gridworld_MDP_Pen_class.east(row,col,Nrow,Ncol,im);
            newrow(3) = row_buf;
            newcol(3) = col_buf;
            newcollision_index(3) = collision_index_buf;
            
            %west
            clear row_buf col_buf collision_index_buf
            [row_buf, col_buf, collision_index_buf]= Gridworld_MDP_Pen_class.west(row,col,Nrow,Ncol,im);
            newrow(4) = row_buf;
            newcol(4) = col_buf;
            newcollision_index(4) = collision_index_buf;
            
        end
    end
    
    methods
        function obj = Gridworld_MDP_Pen_class(ImageFile,eps,targetRow,targetCol,image_type,penalty)
            % construct MDP from image file \ matlab image structure
            % epsilon is transition noise
            if nargin < 5
                image_type = 'file';
            end
            if nargin < 6
                penalty = -10;
            end
            if strcmp(image_type,'file')
                im = imread(ImageFile);
                img = double(rgb2gray(im));
            elseif strcmp(image_type,'matlab')
                img = ImageFile;
            else
                error('unknown image type');
            end
            Nrow = size(img,1);
            Ncol = size(img,2);
            obstacles = find(img == 0);
            nonobstacles = find(img ~= 0);
            
            Ns = Nrow*Ncol;
            Na = 4;
            NumState = Ns+1;
            
            % last state is the die-state
            die_state = NumState;
            
            Pn = zeros(NumState,NumState);  % north
            Ps = zeros(NumState,NumState);  % south
            Pe = zeros(NumState,NumState);  % east
            Pw = zeros(NumState,NumState);  % west
            
            R = -1*ones(NumState,Na);
            
            for i = 1:Nrow
                for j = 1:Ncol
                    %                     current position in state form
                    curpos = sub2ind([Nrow,Ncol],i,j);
                    % generate next positions(north, south, east and west)
                    [newrow,newcol,newcollision_index] = Gridworld_MDP_Pen_class.neigborhood_fn(i,j,Nrow,Ncol,img);
                    %                      next positions in state form
                    for kk = 1:length(newrow)
                        if isfinite(newrow(kk)) == 1 && isfinite(newcol(kk)) == 1
                            % does not hit obstacles
                            nextpos(kk) = sub2ind([Nrow,Ncol],newrow(kk),newcol(kk));
                        else
                            % hit obstacles
                            nextpos(kk) = die_state;
                        end
                    end
                    
                    nextpos_north = nextpos(1);
                    nextpos_south = nextpos(2);
                    nextpos_east = nextpos(3);
                    nextpos_west = nextpos(4);
                    
                    Pn(curpos,nextpos_north) = Pn(curpos,nextpos_north) + 1-eps;
                    Pn(curpos,nextpos_south) = Pn(curpos,nextpos_south) + eps/3;
                    Pn(curpos,nextpos_east) = Pn(curpos,nextpos_east) + eps/3;
                    Pn(curpos,nextpos_west) = Pn(curpos,nextpos_west) + eps/3;
                    
                    Ps(curpos,nextpos_north) = Ps(curpos,nextpos_north) + eps/3;
                    Ps(curpos,nextpos_south) = Ps(curpos,nextpos_south) + 1-eps;
                    Ps(curpos,nextpos_east) = Ps(curpos,nextpos_east) + eps/3;
                    Ps(curpos,nextpos_west) = Ps(curpos,nextpos_west) + eps/3;
                    
                    Pe(curpos,nextpos_north) = Pe(curpos,nextpos_north) + eps/3;
                    Pe(curpos,nextpos_south) = Pe(curpos,nextpos_south) + eps/3;
                    Pe(curpos,nextpos_east) = Pe(curpos,nextpos_east) + 1-eps;
                    Pe(curpos,nextpos_west) = Pe(curpos,nextpos_west) + eps/3;
                    
                    Pw(curpos,nextpos_north) = Pw(curpos,nextpos_north) + eps/3;
                    Pw(curpos,nextpos_south) = Pw(curpos,nextpos_south) + eps/3;
                    Pw(curpos,nextpos_east) = Pw(curpos,nextpos_east) + eps/3;
                    Pw(curpos,nextpos_west) = Pw(curpos,nextpos_west) + 1-eps;
                    
                end
            end
            
            % recurrent state at target
            target = sub2ind([Nrow,Ncol],targetRow,targetCol);
            Pn(target,:) = 0;
            Pn(target,target) = 1;
            Ps(target,:) = 0;
            Ps(target,target) = 1;
            Pe(target,:) = 0;
            Pe(target,target) = 1;
            Pw(target,:) = 0;
            Pw(target,target) = 1;
            
            R(target,:) = 0;
            
            % die_state (state NumState) is recurrent as well
            %              recurrent state at die state
            Pn(die_state,:) = 0;
            Pn(die_state,die_state) = 1;
            Ps(die_state,:) = 0;
            Ps(die_state,die_state) = 1;
            Pe(die_state,:) = 0;
            Pe(die_state,die_state) = 1;
            Pw(die_state,:) = 0;
            Pw(die_state,die_state) = 1;
            
            R(die_state,:) = penalty;
            
            %              add the die_state back to the system
            nonobstacles_total = [nonobstacles; die_state];
            
            % states and rewards only define for nonobstacles
            Pn = Pn(nonobstacles_total,nonobstacles_total);
            Ps = Ps(nonobstacles_total,nonobstacles_total);
            Pe = Pe(nonobstacles_total,nonobstacles_total);
            Pw = Pw(nonobstacles_total,nonobstacles_total);
            
            
            % with four actions: north south east west
            P(:,:,1) = Pn;
            P(:,:,2) = Ps;
            P(:,:,3) = Pe;
            P(:,:,4) = Pw;
            
            % reward for nonobstacles
            R = R(nonobstacles_total,:);
            
            obj@Finite_MDP_class(P,R);
            obj.Nrow = Nrow;
            obj.Ncol = Ncol;
            obj.img = img;
            obj.obstacles = obstacles;
            obj.nonobstacles = nonobstacles;
            obj.targetRow = targetRow;
            obj.targetCol = targetCol;
            obj.eps = eps;
        end
        
        function im = val2image(obj,val)
            im = zeros([obj.Nrow, obj.Ncol]);
            % define the obstacle to be low reward
            % get rid of the fake recurrent state: die_state
            im(obj.nonobstacles) = val;
        end
        
        function [s] = getState(obj,row,col)
            % return state index for row and column
            if row==0 || col==0
                s = numel(obj.nonobstacles)+1; % die state
            else
                s = find(obj.nonobstacles == sub2ind([obj.Nrow,obj.Ncol],row,col),1);
                if isempty(s)
                    s = numel(obj.nonobstacles)+1; % obstacle - go to die state
                end
            end
        end
        
        function [row,col] = getRowCol(obj,s)
            % return row and column for state index
            row = zeros(size(s)); col = row;
            for i = 1: numel(s)
                if s(i) > numel(obj.nonobstacles) % die state
                    row(i) = 0;
                    col(i) = 0;
                else
                    [row(i),col(i)] = ind2sub([obj.Nrow,obj.Ncol],obj.nonobstacles(s(i)));
                end
            end
        end
        
        
        function [row,col] = getPath(obj,pol,row0,col0,maxIter)
            % get path until target starting from (row0,col0)
            pol_im = zeros(obj.Nrow,obj.Ncol);
            pol_im(obj.nonobstacles) = pol;
            
            if nargin < 5
                maxIter = 1e3;
            end
            
            row = zeros(maxIter,1);
            col = zeros(maxIter,1);
            row(1) = row0;
            col(1) = col0;
            for i = 2:maxIter
                if pol_im(row(i-1),col(i-1)) == 1   % north
                    [row(i),col(i)] = Gridworld_MDP_Pen_class.north(row(i-1),col(i-1),obj.Nrow,obj.Ncol,obj.img);
                elseif pol_im(row(i-1),col(i-1)) == 2   % south
                    [row(i),col(i)] = Gridworld_MDP_Pen_class.south(row(i-1),col(i-1),obj.Nrow,obj.Ncol,obj.img);
                elseif pol_im(row(i-1),col(i-1)) == 3   % east
                    [row(i),col(i)] = Gridworld_MDP_Pen_class.east(row(i-1),col(i-1),obj.Nrow,obj.Ncol,obj.img);
                elseif pol_im(row(i-1),col(i-1)) == 4% west
                    [row(i),col(i)] = Gridworld_MDP_Pen_class.west(row(i-1),col(i-1),obj.Nrow,obj.Ncol,obj.img);
                else
                    error('unknown action');
                end
                if isfinite(row(i)) ~= 1 || isfinite(col(i)) ~= 1
                    disp('Hit an obstacle, mission failed')
                    break;
                end
                if obj.targetRow == row(i) && obj.targetCol == col(i)
                    disp('Arrive at target, mission accomplised')
                    break;
                end
            end
            row = row(1:i);
            col = col(1:i);
        end
    end
end