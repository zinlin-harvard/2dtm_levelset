classdef LSGeo
    % controls pmls, periodicity, etc. that are separate from the
    % properties of a LevelSet class. Methods to generate meshes, subdomain
    % keys, boundary keys, etc.
    
    
    properties
        ls; % associated LevelSet object
    end
    
    properties (SetAccess = private)
        bc_xy; % 'bloch', 'per' or 'pec'
        pml_xy; % 0=no, 1=yes
        lpml_xy; % pml lengths
                
        % geometry determined by two rectangles:
        %     level-set rectangle, from p0ls (lower-left corner) to p1ls (upper-right corner)
        %     domain rectangle, from p0 to p1
        p0ls;
        p1ls;
        p0;
        p1;
    end
    
    methods
        function obj = LSGeo(lsi, bc_xyi, pml_xyi, lpml_xyi)
            obj.ls = lsi; % 'shallow' copy (don't delete lsi!)
            obj.bc_xy = bc_xyi;
            obj.pml_xy = pml_xyi;
            obj.lpml_xy = lpml_xyi;
            obj.p0ls = [min(obj.ls.x(:))-0.5*obj.ls.g.dx min(obj.ls.y(:))-0.5*obj.ls.g.dy];
            obj.p1ls = [max(obj.ls.x(:))+0.5*obj.ls.g.dx max(obj.ls.y(:))+0.5*obj.ls.g.dy];
            obj.p0 = obj.p0ls;
            obj.p1 = obj.p1ls;
            for i = 1:2
                if (obj.pml_xy(i))
                    obj.p0(i) = obj.p0(i) - obj.lpml_xy(i);
                    obj.p1(i) = obj.p1(i) + obj.lpml_xy(i);
                end
            end
        end
        
        % edge_length optional, probably should be left at ls.g.dx
        function[] = gen_mesh(obj, mesh_fname, constraint_fname, edge_length)
            if (nargin<4)
                edge_length = min(obj.ls.g.dx, obj.ls.g.dy);
            end
            obj.write_constraint_file(constraint_fname);
            
            % TODO: fix periodicity implementation, right now off (0 0)
            call_str = sprintf('./gen_optimized_mesh %s %s %g %g %g %g %g 0 0', ...
                mesh_fname, constraint_fname, edge_length, ...
                obj.p0(1), obj.p0(2), obj.p1(1), obj.p1(2));
            system(call_str);
            
            call_str = sprintf('./convert_gmsh_to_xml %s mesh.xml', mesh_fname);
            system(call_str);
        end
        
        function[] = write_constraint_file(obj, fname)
            LSGeo.write_isoc_as_constraint_file(obj.ls, fname);
            
            if (obj.pml_xy(1) && obj.pml_xy(2))
                % four sides of inner rectangle
                dlmwrite(fname, [obj.p0ls(1:2) obj.p0ls(1) obj.p1ls(2)], '-append', 'precision', 16, 'delimiter', ' ');
                dlmwrite(fname, [obj.p0ls(1) obj.p1ls(2) obj.p1ls(1:2)], '-append', 'precision', 16, 'delimiter', ' ');
                dlmwrite(fname, [obj.p1ls(1:2) obj.p1ls(1) obj.p0ls(2)], '-append', 'precision', 16, 'delimiter', ' ');
                dlmwrite(fname, [obj.p1ls(1) obj.p0ls(2) obj.p0ls(1:2)], '-append', 'precision', 16, 'delimiter', ' ');
                % connect to corners
                doff = 0; % useful for testing, TODO: remove once periodic-mesh implemented properly
                dlmwrite(fname, [obj.p0(1)+doff obj.p0ls(2) obj.p0ls(1:2)], '-append', 'precision', 16, 'delimiter', ' ');
                dlmwrite(fname, [obj.p0ls(1) obj.p0(2)+doff obj.p0ls(1:2)], '-append', 'precision', 16, 'delimiter', ' ');
                dlmwrite(fname, [obj.p0(1)+doff obj.p1ls(2) obj.p0ls(1) obj.p1ls(2)], '-append', 'precision', 16, 'delimiter', ' ');
                dlmwrite(fname, [obj.p0ls(1) obj.p1(2)-doff obj.p0ls(1) obj.p1ls(2)], '-append', 'precision', 16, 'delimiter', ' ');
                dlmwrite(fname, [obj.p1(1)-doff obj.p1ls(2) obj.p1ls(1:2)], '-append', 'precision', 16, 'delimiter', ' ');
                dlmwrite(fname, [obj.p1ls(1) obj.p1(2)-doff obj.p1ls(1:2)], '-append', 'precision', 16, 'delimiter', ' ');
                dlmwrite(fname, [obj.p1(1)-doff obj.p0ls(2) obj.p1ls(1) obj.p0ls(2)], '-append', 'precision', 16, 'delimiter', ' ');
                dlmwrite(fname, [obj.p1ls(1) obj.p0(2)+doff obj.p1ls(1) obj.p0ls(2)], '-append', 'precision', 16, 'delimiter', ' ');
            elseif (obj.pml_xy(1))
                dlmwrite(fname, [obj.p0ls(1:2) obj.p0ls(1) obj.p1ls(2)], '-append', 'precision', 16, 'delimiter', ' ');
                dlmwrite(fname, [obj.p1ls(1:2) obj.p1ls(1) obj.p0ls(2)], '-append', 'precision', 16, 'delimiter', ' ');
            elseif (obj.pml_xy(2))
                dlmwrite(fname, [obj.p0ls(1) obj.p1ls(2) obj.p1ls(1:2)], '-append', 'precision', 16, 'delimiter', ' ');
                dlmwrite(fname, [obj.p1ls(1) obj.p0ls(2) obj.p0ls(1:2)], '-append', 'precision', 16, 'delimiter', ' ');
            end
        end
        
        % labels: 1 = interior of level-set object
        %         0 = exterior of level-set object (embedding)
        %         2 = pml_x
        %         3 = pml_y
        %         5 = pml_xy
        function[domains] = mark_domains_xml(obj, mesh, fname)
            x_mesh = mesh.centroids(:,1);
            y_mesh = mesh.centroids(:,2);
            phi_mesh = obj.ls.ls_interp(x_mesh, y_mesh);
            domains = zeros(size(phi_mesh));
            domains = domains + 2.*(x_mesh<obj.p0ls(1)) + 2.*(x_mesh>obj.p1ls(1)) ...
                              + 3.*(y_mesh<obj.p0ls(2)) + 3.*(y_mesh>obj.p1ls(2));
            domains(phi_mesh<0) = 1;
            mesh.viewPatch(domains);
            
            fid = fopen(fname,'w');
			fprintf(fid, '<?xml version="1.0"?>\n');
            fprintf(fid, '<dolfin xmlns:dolfin="http://fenicsproject.org">\n');
            fprintf(fid, '  <mesh_function>\n');
            fprintf(fid, '    <mesh_value_collection name="f" type="uint" dim="2" size="%i">\n', mesh.numTri);
            for i = 1:mesh.numTri
                fprintf(fid, '      <value cell_index="%i" local_entity="0" value="%i" />\n', i-1, domains(i));
            end
            fprintf(fid, '    </mesh_value_collection>\n');
            fprintf(fid, '  </mesh_function>\n');
            fprintf(fid, '</dolfin>\n');
			fclose(fid);
        end
        
        % use domains argument if want material boundaries to be identified
        function[facet_bnds] = mark_bnds_xml(obj, mesh, fname, domains)
            tr = triangulation(mesh.t, mesh.v);
            
            % most facets have continuity bc
            facet_bnds = ones(mesh.numTri, 3);
            
            % free-boundary edges
            edges_fb = freeBoundary(tr);
            
            % labels: 0 = pec, 1 = interior, 2 = bloch_x (= per_x), 3 =
            % bloch_y (= per_y), 4 = material boundary (if requested)
            vi = [2 3; 1 3; 1 2]; % labeling of edges found by trial-and-error
            t_nz = cell2mat(edgeAttachments(tr, edges_fb));
            for ti = t_nz.'
                for i = 1:3
                    edge_i = mesh.t(ti, vi(i,:));
                    if (ismember(edge_i, edges_fb))
                        x_dist = norm(mesh.v(edge_i(1),1) - mesh.v(edge_i(2),1));
                        y_dist = norm(mesh.v(edge_i(1),2) - mesh.v(edge_i(2),2));
                        if (x_dist<y_dist && (strcmp(obj.bc_xy(1), 'bloch') || strcmp(obj.bc_xy(1), 'per')))
                            facet_bnds(ti,i) = 2;
                        elseif (y_dist<x_dist && (strcmp(obj.bc_xy(2), 'bloch') || strcmp(obj.bc_xy(2), 'per')))
                            facet_bnds(ti,i) = 3;
                        else
                            facet_bnds(ti,i) = 0;
                        end
                    end
                end
            end
            
            % material boundaries, if desired
            if (nargin==4)
                %  need to find:
                %    all triangles with domain=1 and neighbors domain=0, and
                %    all triangles with domain=0 and neighbors domain=1
                %    (because same edge belongs to two triangles)
                tri = (1:mesh.numTri).'; 
                for i = [0 1]
                    interior_tri = tri(domains==1-i);
                    nbr_tri = tr.neighbors(interior_tri);
                    bnd_tri1 = interior_tri(domains(nbr_tri(:,1))==i);
                    bnd_tri2 = interior_tri(domains(nbr_tri(:,2))==i);
                    bnd_tri3 = interior_tri(domains(nbr_tri(:,3))==i);
                    facet_bnds(bnd_tri1,1) = 4;
                    facet_bnds(bnd_tri2,2) = 4;
                    facet_bnds(bnd_tri3,3) = 4;
                end
            end
            
            fid = fopen(fname,'w');
			fprintf(fid, '<?xml version="1.0"?>\n');
            fprintf(fid, '<dolfin xmlns:dolfin="http://fenicsproject.org">\n');
            fprintf(fid, '  <mesh_function>\n');
            fprintf(fid, '    <mesh_value_collection name="f" type="uint" dim="1" size="%i">\n', 3*mesh.numTri);
            for i = 1:mesh.numTri
                for j = 1:3
                    fprintf(fid, '      <value cell_index="%i" local_entity="%i" value="%i" />\n', i-1, j-1, facet_bnds(i,j));
                end
            end
            fprintf(fid, '    </mesh_value_collection>\n');
            fprintf(fid, '  </mesh_function>\n');
            fprintf(fid, '</dolfin>\n');
			fclose(fid);
        end
    end
    
    methods (Static)
        
        function[] = write_isoc_as_constraint_file(ls, constraint_fname)
            c = (contourc(ls.phi, [0 0]).' - 1) * ls.g.dx + ls.g.x0;
            
            % first row of c is: contour (i.e. 0), num_points (last is redundant)
            c_lines = [c(2:end-1,:) c(3:end,:)];
            
            % remove very short lines
            sl_tol = 1e-6;
            line_lengths = sqrt((c_lines(:,3) - c_lines(:,1)).^2 + (c_lines(:,4) - c_lines(:,2)).^2);
            cond = line_lengths<sl_tol;
            while (sum(cond)>0)
                sl_ind = find(cond, 1); % index of first short line
                N = size(c_lines,1);
                prev_ind = (sl_ind-1) * (sl_ind>1) + N * (sl_ind==1);
                next_ind = (sl_ind+1) * (sl_ind<N) + 1 * (sl_ind==N);
                prev_point = c_lines(prev_ind, 3:4);
                next_point = c_lines(next_ind, 1:2);
                new_point = 0.5 * (prev_point + next_point);
                c_lines(prev_ind, 3:4) = new_point;
                c_lines(next_ind, 1:2) = new_point;
                c_lines(sl_ind, :) = [];
                line_lengths = sqrt((c_lines(:,3) - c_lines(:,1)).^2 + (c_lines(:,4) - c_lines(:,2)).^2);
                cond = line_lengths<sl_tol;
                fprintf('deleted short line\n');
            end
            
            % write to file
            dlmwrite(constraint_fname, c_lines, 'precision', 16, 'delimiter', ' ');
        end   
    end
    
end

