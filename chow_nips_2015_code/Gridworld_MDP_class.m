classdef Gridworld_MDP_class < Finite_MDP_class
    % Gridworld domain with obstacles. Actions are to
    % {north,south,east,west} and w.p. epsilon a random movement is
    % performed (noisy transitions).
    properties
        Nrow = 1;      % image rows
        Ncol = 1;      % image columns
        img = [];      % image
        obstacles = [];% indices of obtacles in image
        non_obstacles; % indices of obtacles in image
        targetRow = 1;
        targetCol = 1;
    end
    methods (Static)
        function [newrow,newcol] = north(row,col,Nrow,Ncol,im)
            newrow = max(row-1,1);
            newcol = col;
            if im(newrow,newcol) == 0   % obstacle
                newrow = row;
                newcol = col;
            end
        end
        function [newrow,newcol] = south(row,col,Nrow,Ncol,im)
            newrow = min(row+1,Nrow);
            newcol = col;
            if im(newrow,newcol) == 0   % obstacle
                newrow = row;
                newcol = col;
            end
        end
        function [newrow,newcol] = east(row,col,Nrow,Ncol,im)
            newrow = row;
            newcol = min(col+1,Ncol);
            if im(newrow,newcol) == 0   % obstacle
                newrow = row;
                newcol = col;
            end
        end
        function [newrow,newcol] = west(row,col,Nrow,Ncol,im)
            newrow = row;
            newcol = max(col-1,1);
            if im(newrow,newcol) == 0   % obstacle
                newrow = row;
                newcol = col;
            end
        end
        function [rows,cols] = neighbors(row,col,Nrow,Ncol,im)
            [rows,cols] = Gridworld_MDP_class.north(row,col,Nrow,Ncol,im);
            [newrow,newcol] = Gridworld_MDP_class.south(row,col,Nrow,Ncol,im);
            rows = [rows,newrow]; cols = [cols,newcol];
            [newrow,newcol] = Gridworld_MDP_class.east(row,col,Nrow,Ncol,im);
            rows = [rows,newrow]; cols = [cols,newcol];
            [newrow,newcol] = Gridworld_MDP_class.west(row,col,Nrow,Ncol,im);
            rows = [rows,newrow]; cols = [cols,newcol];
        end
    end
    methods 
        function obj = Gridworld_MDP_class(ImageFile,eps,targetRow,targetCol)
            % construct MDP from image file
            % epsilon is transition noise
            im = imread(ImageFile);
            img = double(rgb2gray(im));
            Nrow = size(img,1);
            Ncol = size(img,2);
            obstacles = find(img == 0);
            non_obstacles = find(img ~= 0);
            target = sub2ind([Nrow,Ncol],targetRow,targetCol);
            Ns = Nrow*Ncol;
            Na = 4;
            Pn = zeros(Ns,Ns);  % north
            Ps = zeros(Ns,Ns);  % south
            Pe = zeros(Ns,Ns);  % east
            Pw = zeros(Ns,Ns);  % west
            R = -1*ones(Ns,Na);
            R(target,:) = 0;
            for row = 1:Nrow
                for col = 1:Ncol
                    curpos = sub2ind([Nrow,Ncol],row,col);
                    [rows,cols] = Gridworld_MDP_class.neighbors(row,col,Nrow,Ncol,img);
                    neighbor_inds = sub2ind([Nrow,Ncol],rows,cols);
                    Pn(curpos,neighbor_inds(1)) = Pn(curpos,neighbor_inds(1)) + 1-eps;
                    Pn(curpos,neighbor_inds(2)) = Pn(curpos,neighbor_inds(2)) + eps/3;
                    Pn(curpos,neighbor_inds(3)) = Pn(curpos,neighbor_inds(3)) + eps/3;
                    Pn(curpos,neighbor_inds(4)) = Pn(curpos,neighbor_inds(4)) + eps/3;
                    Ps(curpos,neighbor_inds(2)) = Ps(curpos,neighbor_inds(2)) + 1-eps;
                    Ps(curpos,neighbor_inds(1)) = Ps(curpos,neighbor_inds(1)) + eps/3;
                    Ps(curpos,neighbor_inds(3)) = Ps(curpos,neighbor_inds(3)) + eps/3;
                    Ps(curpos,neighbor_inds(4)) = Ps(curpos,neighbor_inds(4)) + eps/3;
                    Pe(curpos,neighbor_inds(3)) = Pe(curpos,neighbor_inds(3)) + 1-eps;
                    Pe(curpos,neighbor_inds(1)) = Pe(curpos,neighbor_inds(1)) + eps/3;
                    Pe(curpos,neighbor_inds(2)) = Pe(curpos,neighbor_inds(2)) + eps/3;
                    Pe(curpos,neighbor_inds(4)) = Pe(curpos,neighbor_inds(4)) + eps/3;
                    Pw(curpos,neighbor_inds(4)) = Pw(curpos,neighbor_inds(4)) + 1-eps;
                    Pw(curpos,neighbor_inds(1)) = Pw(curpos,neighbor_inds(1)) + eps/3;
                    Pw(curpos,neighbor_inds(2)) = Pw(curpos,neighbor_inds(2)) + eps/3;
                    Pw(curpos,neighbor_inds(3)) = Pw(curpos,neighbor_inds(3)) + eps/3;
                end
            end
            Pn(target,:) = 0*Pn(target,:); Pn(target,target) = 1;
            Ps(target,:) = 0*Ps(target,:); Ps(target,target) = 1;
            Pe(target,:) = 0*Pe(target,:); Pe(target,target) = 1;
            Pw(target,:) = 0*Pw(target,:); Pw(target,target) = 1;
            Pn = Pn(non_obstacles,:); Pn = Pn(:,non_obstacles);
            Ps = Ps(non_obstacles,:); Ps = Ps(:,non_obstacles);
            Pe = Pe(non_obstacles,:); Pe = Pe(:,non_obstacles);
            Pw = Pw(non_obstacles,:); Pw = Pw(:,non_obstacles);
            R = R(non_obstacles,:);
            P = cat(3,Pn,Ps,Pe,Pw);
            obj@Finite_MDP_class(P,R);
            obj.Nrow = Nrow;
            obj.Ncol = Ncol;
            obj.img = img;
            obj.obstacles = obstacles;
            obj.non_obstacles = non_obstacles;
            obj.targetRow = targetRow;
            obj.targetCol = targetCol;
        end
        function [im] = val2image(obj,val)
            % put values (for states) on the image
            im = zeros(obj.Nrow,obj.Ncol);
            im(obj.non_obstacles) = val;
        end
        function [row,col] = getPath(obj,pol,row0,col0,maxIter)
            % get path until target starting from (row0,col0)
            pol_im = zeros(obj.Nrow,obj.Ncol);
            pol_im(obj.non_obstacles) = pol;
            
             if nargin < 5
                 maxIter = 1e3;
             end
            
            row = zeros(maxIter,1);
            col = zeros(maxIter,1);
            row(1) = row0;
            col(1) = col0;
            for i = 2:maxIter
                if pol_im(row(i-1),col(i-1)) == 1   % north
                    [row(i),col(i)] = Gridworld_MDP_class.north(row(i-1),col(i-1),obj.Nrow,obj.Ncol,obj.img);
                elseif pol_im(row(i-1),col(i-1)) == 2   % south
                    [row(i),col(i)] = Gridworld_MDP_class.south(row(i-1),col(i-1),obj.Nrow,obj.Ncol,obj.img);
                elseif pol_im(row(i-1),col(i-1)) == 3   % east
                    [row(i),col(i)] = Gridworld_MDP_class.east(row(i-1),col(i-1),obj.Nrow,obj.Ncol,obj.img);
                elseif pol_im(row(i-1),col(i-1)) == 4% west
                    [row(i),col(i)] = Gridworld_MDP_class.west(row(i-1),col(i-1),obj.Nrow,obj.Ncol,obj.img);
                else
                    error('unknown action');
                end
                if obj.targetRow == row(i) && obj.targetCol == col(i) 
                    break;
                end
            end
            row = row(1:i);
            col = col(1:i);
        end
    end
end