classdef Mesh
	% Mesh: for a given triangulation, determine important information
	% (i.e. obj.volume, surface area, normals), and calculate derivatives, etc.
	% TODO: TriRep --> triangulation
	properties
		
		% vertices/triangles
		v;
		t;
		numVert;
		numTri;
		vNbrs; % N x m cell array containing neighbors for each vertex
		
		% characteristic quantities
		vol;
		SA;
		
		% normals
		vertNorms;
		triNorms;
		centroids; % for viewing triNorms
		panelAreas;
		
		% mesh quality
		eta;
		gamma;
		etaM; % mean(eta)
		gammaM; % mean(gamma)
		
		% verbose?
		verbose = 1;
		
	end
	
	methods
		% Constructor
		%	Accepts either: 
		%		one argument, filename (must end in "off" or "gmsh")
		%		or two arguments, vertices + triangles
		function obj = Mesh(varargin)
			offType = 0;
			if(length(varargin)==1)
				[~,~,ext] = fileparts(varargin{1});
				if(strcmp(ext,'.msh'))
					[obj.v, obj.t] = Mesh.LoadGmshFile(varargin{1});
				elseif(strcmp(ext,'.off'))
					[obj.v, obj.t] = Mesh.LoadOffFile(varargin{1});
					offType = 1;
				else
					error('Incorrect filetype extension: %s', ext);
				end
			else
				obj.v = varargin{1};
				obj.t = varargin{2};
			end
			obj.numVert = size(obj.v,1);
			obj.numTri = size(obj.t,1);
			obj = obj.processMesh;
			
			% correct if normals facing inward
			if(offType)
				if(obj.vol<0)
					obj.t = [obj.t(:,1) obj.t(:,3) obj.t(:,2)];
					warning('normals facing wrong direction in OFF file');
					obj = obj.processMesh;
				end
			end
		end
		
		% Copy constructor
		function obj = MeshCopy(m0)
			obj = Mesh(m0.v, m0.t);
		end
		
		% Set verbosity
		function obj = setVerbose(obj, v)
			obj.verbose = (v>=1); % 1 if v>=1, else 0
		end
		
		% move vertices and process mesh at same time
		function obj = moveV(obj, dv)
			if( (size(dv,1)~=obj.numVert) )
				error('dv has incorrect number of vertices\n');
			end
			if( size(dv,2)==1 )
				% multiply by vertNorms
				dv = repmat(dv, [1 3]).*obj.vertNorms;
			elseif( size(dv,2)~=3 )
				error('dv has incorrect number of columns\n');
			end
			obj.v = obj.v + dv;
			obj = obj.processMesh;
		end
		
		% Compute vol, surface area, and vertex and panel normal vectors
		function[obj] = processMesh(obj)
			obj.vol = 0;
			obj.SA = 0;
			obj.triNorms = zeros(obj.numTri, 3);			
			obj.centroids = zeros(obj.numTri, 3);
			obj.panelAreas = zeros(obj.numTri, 1);
			obj.vertNorms = zeros(obj.numVert, 3);
            tic;
            i = 1:obj.numTri;
            x1 = obj.v( obj.t(i,1),:);
            x2 = obj.v( obj.t(i,2),:);
            x3 = obj.v( obj.t(i,3),:);
            a = x2-x1;
            b = x3-x1;
            c = x3-x2;
            cross12 = cross(a,b);
            norm_a = vnorm2(a);
            norm_b = vnorm2(b);
            norm_c = vnorm2(c);
            norm_cross12 = sqrt(sum(cross12.^2, 2));
            
            obj.vol = (1/6) * sum( dot(x1, cross(x2, x3), 2) );
            obj.SA = 0.5 * sum( norm_cross12 );
            
            obj.centroids(i,:) = (1/3) * (x1+x2+x3);
            obj.triNorms(i,:) = cross12 ./ repmat(norm_cross12, 1, 3);
            obj.panelAreas(i) = 0.5 * norm_cross12;
            
            % angles
            a1 = repmat(acos( dot(a,b,2) ./ (norm_a.*norm_b) ), 1, 3);
            a2 = repmat(acos( dot(-a,c,2) ./ (norm_a.*norm_c) ), 1, 3);
            a3 = repmat(acos( dot(-b,-c,2) ./ (norm_b.*norm_c) ), 1, 3);
            
            obj.vertNorms( obj.t(i,1),:) = a1.*obj.triNorms(i,:);
            obj.vertNorms( obj.t(i,2),:) = a2.*obj.triNorms(i,:);
            obj.vertNorms( obj.t(i,3),:) = a3.*obj.triNorms(i,:);
            
			obj.vertNorms = real(obj.vertNorms);
			obj.vertNorms = obj.vertNorms./repmat(vnorm2(obj.vertNorms), [1 3]);
			
            obj = obj.meshQuality();
            
			if(obj.verbose)
				fprintf('  Volume: %g \n', obj.vol);
				fprintf('  Surface Area: %g \n', obj.SA);
				fprintf('  Eta: %g\n  Gamma: %g\n', obj.etaM, obj.gammaM);  
			end
        end
        
        function[obj] = reorder_triangles(obj)
            for i = 1:obj.numTri
                obj.t(i,:) = sort(obj.t(i,:), 'ascend');
            end
        end
        
		% Set vNbrs cell array
		function[obj] = setNbrs(obj, tr)
			tNbrs = vertexAttachments(tr);
			obj.vNbrs = cell(obj.numVert,1);
			for i=1:obj.numVert
				vn = obj.t(tNbrs{i}.',:);
				obj.vNbrs{i} = setdiff(unique(vn(:)), i);
			end	
		end
		
		% measure mesh quality (using same eta, gamma as defined in GMSH)
		function[obj] = meshQuality(obj)
			a = zeros(obj.numTri,3);
			v1 = obj.v(obj.t(:,1),:);
			v2 = obj.v(obj.t(:,2),:);
			v3 = obj.v(obj.t(:,3),:);
			sin1 = vnorm2(cross(v2-v1,v3-v1,2)) ./ vnorm2(v2-v1) ./ vnorm2(v3-v1);
			sin2 = vnorm2(cross(v3-v2,v1-v2,2)) ./ vnorm2(v3-v2) ./ vnorm2(v1-v2);
			sin3 = vnorm2(cross(v1-v3,v2-v3,2)) ./ vnorm2(v1-v3) ./ vnorm2(v2-v3);
			a(:,1) = 180/pi*atan2( vnorm2(cross(v2-v1,v3-v1,2)), dot(v2-v1,v3-v1,2) );
			a(:,2) = 180/pi*atan2( vnorm2(cross(v3-v2,v1-v2,2)), dot(v3-v2,v1-v2,2) );
			a(:,3) = 180/pi*atan2( vnorm2(cross(v1-v3,v2-v3,2)), dot(v1-v3,v2-v3,2) );
			obj.eta = min(a/60,[],2);
			obj.gamma = 4.*sin1.*sin2.*sin3./(sin1+sin2+sin3);
			obj.etaM = mean(obj.eta);
			obj.gammaM = mean(obj.gamma);
		end
		
		% derivative of volume with respect to each vertex, in the normal
		%	direction
		function [dV] = dVdxn(obj)
			dV = zeros(obj.numVert, 1);
			% Loop through triangles, add new deriv term to previous ones
			for i = 1:obj.numTri
				% extract triangle info
				v1 = obj.t(i,1);
				v2 = obj.t(i,2);
				v3 = obj.t(i,3);
				x1 = obj.v(v1,:);
				x2 = obj.v(v2,:);
				x3 = obj.v(v3,:);
				
				% vertices
				dV(v1) = dV(v1) + (1/6)*dot( obj.vertNorms(v1,:), cross(x2,x3) );
				dV(v2) = dV(v2) + (1/6)*dot( obj.vertNorms(v2,:), cross(x3,x1) );
				dV(v3) = dV(v3) + (1/6)*dot( obj.vertNorms(v3,:), cross(x1,x2) );
			end
		end
		
		% derivative of surface area with respect to each vertex
		function[dA] = dAdxn(obj)
			dA = zeros(obj.numVert, 1);
			for i=1:obj.numTri
				v1 = obj.t(i,1);
				v2 = obj.t(i,2);
				v3 = obj.t(i,3);
				x1 = obj.v(v1,:);
				x2 = obj.v(v2,:);
				x3 = obj.v(v3,:);
				
				% deriv for v1
				a = x1 - x3; b = x2 - x3;
				nab = cross(a,b)/norm(cross(a,b));
				dA(v1) = dA(v1) + 0.5*dot( cross(obj.vertNorms(v1,:), b), nab);
				
				% deriv for v1
				a = x2 - x1; b = x3 - x1;
				nab = cross(a,b)/norm(cross(a,b));
				dA(v1) = dA(v1) + 0.5*dot( cross(obj.vertNorms(v2,:), b), nab);
				
				% deriv for v1
				a = x3 - x1; b = x2 - x1;
				nab = cross(a,b)/norm(cross(a,b));
				dA(v1) = dA(v1) + 0.5*dot( cross(obj.vertNorms(v3,:), b), nab);
			end
		end
		
		% Calculate total angular defect (=4pi if homeomorphic to sphere)
		% angDefect contains angular defect at each vertex
		function [angDefect, defectVec] = calcDefect(obj)
			defectVec = 2*pi*ones(obj.numVert,1);
			for i=1:obj.numTri
				tri = obj.t(i,:);
				v1 = tri(1); v2 = tri(2); v3 = tri(3);
				x1 = obj.v(v1,:); x2 = obj.v(v2,:); x3 = obj.v(v3,:);
				a = x2 - x1; b = x3 - x1; c = x3 - x2;
				
				defectVec(tri(1)) = defectVec(tri(1)) - acos( dot(a,b) / norm(a) / norm(b) );
				defectVec(tri(2)) = defectVec(tri(2)) - acos( dot(-a,c) / norm(a) / norm(c) );
				defectVec(tri(3)) = defectVec(tri(3)) - acos( dot(b,c) / norm(b) / norm(c) );
			end
			angDefect = sum(defectVec);
		end
		
		% Test topology (i.e. no edge shared by 3 panels)
		function topologyOK = testTopology(obj)
			topologyOK = true;
			tr = TriRep(obj.t, obj.v);
			e = edges(tr);
			for i=1:size(e,1)
				triE = edgeAttachments(tr, e(i,:));
				if(length(triE{1})~=2)
					topologyOK = false;
					break;
				end
			end
		end

		% calculate all distances between neighbors and non-neighbors
		function[dn, dnn, df] = nbrDist(obj, delta, plotYN)
			if(nargin==1), delta = 1; end
			if(nargin<3), plotYN = 1; end;
			dn = [];
			dnn = [];
			df = [];
            tr = TriRep(obj.t, obj.v);
            obj = obj.setNbrs(tr);

			for i=1:obj.numVert
				% find all neighbors
				nbrVerts = obj.vNbrs{i};
				
				% find the next-neighbors
				nnbrVerts = [];
				for j=1:length(nbrVerts)
					nnbrVerts = [nnbrVerts; obj.vNbrs{nbrVerts(j)}];
				end
				nnbrVerts = unique(nnbrVerts);
				nnbrVerts = setdiff(nnbrVerts, union(nbrVerts, i));
				
				nbrVerts = nbrVerts(nbrVerts>i);
				for j=1:length(nbrVerts)
					dn = [dn; vnorm2(obj.v(i,:)- obj.v(nbrVerts(j),:))];
				end
				
				
				nnbrVerts = nnbrVerts(nnbrVerts>i);
				for j=1:length(nnbrVerts)
					dnn = [dnn; vnorm2(obj.v(i,:)- obj.v(nnbrVerts(j),:))];
				end
				
				allVerts = (i+1:obj.numVert);
				fVerts = setdiff(allVerts, union(nbrVerts,nnbrVerts));
				for j=1:length(fVerts)
					df = [df; vnorm2(obj.v(i,:)-obj.v(fVerts(j),:))];
				end
			end
			
			fprintf('min nearest neighbor distance: %g\n',min(dn)/delta);
			fprintf('max nearest neighbor distance: %g\n',max(dn)/delta);
			fprintf('min next-nearest neighbor distance: %g\n',min(dnn)/delta);
			fprintf('max next-nearest neighbor distance: %g\n',max(dnn)/delta);
			fprintf('min non-neighbor distance: %g\n',min(df)/delta);
			fprintf('max non-neighbor distance: %g\n',max(df)/delta);

			% plot neighbors and non-neighbors
			if(plotYN)
				dn2 = hist(dn,25);
				dnn2 = hist(dnn,25);
				df2 = hist(df,25);
				figure;
				hist(repmat(dn, [floor(max(df2)/max(dn2)) 1]),25);
				hold on; 
				hist(repmat(dnn, [floor(max(df2)/max(dnn2)) 1]),25);
				hist(df,25);
				h = findobj(gca,'Type','patch');
				set(h(1), 'FaceColor','r'); 
				set(h(2),'FaceColor','g');
				set(h(3),'FaceColor','b');
				legend('neighbors', 'next-neighbors', 'non-neighbors');
				set(gcf,'color','white');
				set(gca,'FontSize',16);
				str = ['nbr mean = ',num2str(mean(dn)), ...
					' nnbr mean = ',num2str(mean(dnn)), ...
					' far mean = ',num2str(mean(df))];
% 				annotation('textbox',[0.15 0.8 0.38 0.1],'string',str,'FontSize',16);
				annotation('textbox',[0.5 0.125 0.38 0.15],'string',str,'FontSize',16,'Color','w');
			end
		end
		
		% Calculate total area of all triangles containing each vertex
		function[areas] = calcVertAreas(obj)
			areas = zeros(obj.numVert,1);
			for i=1:obj.numTri
				v1 = obj.t(i,1);
				v2 = obj.t(i,2);
				v3 = obj.t(i,3);
				x1 = obj.v(v1,:);
				x2 = obj.v(v2,:);
				x3 = obj.v(v3,:);
				
				triArea = 0.5 * norm(cross(x2-x1,x3-x1));			
				areas(v1) = areas(v1) + triArea;
				areas(v2) = areas(v2) + triArea;
				areas(v3) = areas(v3) + triArea;
			end
		end
		
		% view vertex normals
		function[] = viewVertNorms(obj, q)
			if(nargin==1)
				q = obj.vertNorms;
			elseif(size(q,2)==1)
				q = repmat(q,[1 3]).*obj.vertNorms;
			end
			ind = 1:obj.numVert;
			posInd = ind(dot(q,obj.vertNorms,2)>=0);
			negInd = ind(dot(q,obj.vertNorms,2)<0);
			obj.viewMesh;
			hold on;
			quiver3(obj.v(posInd,1), obj.v(posInd,2), obj.v(posInd,3), ...
				q(posInd,1), q(posInd,2), q(posInd,3));
			quiver3(obj.v(negInd,1), obj.v(negInd,2), obj.v(negInd,3), ...
				-q(negInd,1), -q(negInd,2), -q(negInd,3), 'k');
		end
		
		% view triangle normals
		function[] = viewTriNorms(obj, q)
			if(nargin==1)
				q = obj.triNorms;
			elseif(size(q,2)==1)
				q = repmat(q,[1 3]).*obj.triNorms;
			end
			obj.viewMesh;
			hold on;
			quiver3(obj.centroids(:,1), obj.centroids(:,2), obj.centroids(:,3), ...
				q(:,1), q(:,2), q(:,3));
		end
		
		function[p] = viewPatch(obj, q)
			p = patch('Faces',obj.t,'Vertices',obj.v,'FaceVertexCData',q,...
				'FaceColor','flat','CDataMapping','scaled','EdgeColor','none');
			view(3); axis equal; 
			set(gcf,'color','white');
		end
		
		
		% view mesh (colored edges)
		function[] = viewTriMesh(obj)
			trimesh( TriRep(obj.t, obj.v) );
		end
		
		% view mesh (paneled)
		function[] = viewMesh(obj)
			fvs = struct();
			fvs.vertices = obj.v;
			fvs.faces = obj.t;

			p = patch(fvs);
			set(p,'FaceColor',[1 0.25 0]);
% 			set(p,'FaceColor',[1 1 1]*0.7);
			daspect([1,1,1])
			view(3); axis equal
			camlight
			lighting none
			set(gcf,'color','white');
		end
		
		% view derivatives
		function[] = viewDerivs(obj, grad)
			obj.viewMesh;
			hold on; 
			gpos = grad(grad>=0);
			gneg = grad(grad<0);
			vnpos = obj.vertNorms(grad>=0,:);
			vnneg = obj.vertNorms(grad<0,:);
			gpos = repmat(gpos, [1 3]).*vnpos;
			gneg = -repmat(gneg, [1 3]).*vnneg; % flip sign to point outward
			vpos = obj.v(grad>=0,:);
			vneg = obj.v(grad<0,:);
			quiver3(vpos(:,1), vpos(:,2), vpos(:,3), ...
				gpos(:,1), gpos(:,2), gpos(:,3),'b');
			quiver3(vneg(:,1), vneg(:,2), vneg(:,3), ...
				gneg(:,1), gneg(:,2), gneg(:,3),'k');
		end
		
		% Write triangulation data to a GMSH file
		function[] = writeGmshFile(obj, filename)
			fid = fopen(filename,'w'); % Deletes previous contents
			fprintf(fid, '$MeshFormat \n2.2   0   8 \n$EndMeshFormat \n$Nodes \n%g \n', obj.numVert);
			for i = 1:obj.numVert
				fprintf(fid, '%i %4.16g %4.16g %4.16g \n', i, obj.v(i,1), obj.v(i,2), obj.v(i,3));
			end
			fprintf(fid, '$EndNodes \n$Elements \n%i \n', obj.numTri);
			for i = 1:obj.numTri
				fprintf(fid, '%i %i %i %i %i %4.16g %4.16g %4.16g \n', ...
					i, 2, 2, 2, 2, obj.t(i,1), obj.t(i,2), obj.t(i,3));
			end
			fprintf(fid, '$EndElements \n');
			fclose(fid);	
		end
	end
	
	methods (Static)
		
		% check if can delete edge from triangle?
		function[del] = canDelete(edge, tr)
			del = 1;
			aTris = vertexAttachments(tr, edge(1)); aTris = aTris{1};
			bTris = vertexAttachments(tr, edge(2)); bTris = bTris{1};
			for i=1:size(aTris,1)
				for j=1:size(bTris,1)
					if( isempty( setxor(aTris(i,:),bTris(j,:)) ) )
						del = 0;
						break;
					end
				end
			end
		end
		
		% Load Gmsh file into vertices v and triangulation t
		function[v, t] = LoadGmshFile(filename)			
			% open file for reading
			fid = fopen(filename,'r');
			
			% get number of vertices
			notFound = 1;
			numV = [];
			while( notFound && ~feof(fid) )
				lineIn = fgets(fid);
				if( strcmp(lineIn(1:min(6,length(lineIn))),'$Nodes') )
					notFound = 0;
					numV = str2double(fgets(fid));
				end
			end
			if(isempty(numV))
				error('Reached end of file without finding the number of nodes');
			end
			
			% get vertex data
			v = fscanf(fid, '%i %g %g %g', [4 numV]);
			v = v(2:4,:).';
			
			% get number of triangles
			notFound = 1;
			numT = [];
			while( notFound && ~feof(fid) )
				lineIn = fgets(fid);
				if( strcmp(lineIn(1:min(9,length(lineIn))),'$Elements') )
					notFound = 0;
					numT = str2double(fgets(fid));
				end
			end
			if(isempty(numT))
				error('Reached end of file without finding the number of triangles');
			end
			
			% get triangulation data
			t = fscanf(fid, '%i %i %i %i %i %i %i %i', [8 numT]);
			t = t(6:8,:).';
			
			fclose(fid);
		end
		
		% import an OFF file as would be the output from a CGAL command
		function [v,t] = LoadOffFile(filename)
			% open file for reading
			fid = fopen(filename,'r');
			
			% get number of vertices and triangles
			fgets(fid); % skip the first line ("OFF")
			data = fscanf(fid, '%i %i %i', [3 1]);
			numV = data(1);
			numT = data(2);
			
			% get vertices
			fgets(fid); % skip a blank line
			v = fscanf(fid, '%g %g %g', [3 numV]).';
			
			% get triangles
			t = fscanf(fid, '%i %i %i %i', [4 numT]);
			t = t([2,3,4],:).' + 1; % OFF uses 0-based indexing, inverted definition of normals
			
			fclose(fid);
		end
		
	end
end

function[val] = vnorm2(v)
    val = sqrt(sum(v.^2,2));
end