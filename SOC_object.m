classdef SOC_object
    properties
        %parameters of an object
        
        Rho = 1/0.28*1.68e-8 ;%Resistivity of the object in Ohm meter
        H = 1e-3;%Height of the object in meter
        D1 = 1e-3;%Inner diameter of the object in meter
        D2 = 1.5e-3;%Outer diameter of the object in meter
        xpos = 1e-3;%Location of the object in x in meter
        zpos = 1e-3;%Location of the object in z in meter
        
        %parameters specific for a coil
        
        N = 100;%Number of windings of the coil
        is_coil = 0;%Boolean that specifies that this is a coil
        is_circular = 0;%Boolean that specifies that this is a toroidal coil
        
        %parameters specific for a nozzle
        
        taper_width = 0;%Horizontal distance of the taper of the nozzle in meter
        taper_angle = 0;%Angle of the taper of the nozzle in radians
        hole_offset = 0;%How much the hole of the nozzle is offset from the center in meter
        rotation = 1;%Rotation of the nozzle in degrees
        
        %"mesh" parameters
        
        nx = 11;%Number of elements on the bottom/top boundary mesh in the x direction
        nz = 10;%Number of elements on the left/right/center boundary mesh in the z direction
        nphi = 20;%Number of elements of the coil in radial direction
        nxb = 23;%Total number of elements on the left and right side in the x direction
        nxbl = 3;%Number of elements on the left side in the x direction
        nxbr = 20;%Number of elements on the right side in the x direction
        nxc = 3;%Number of elements on the center in the x direction
        nzb = 21;%Total number of elements on the top and bottom side in the z direction
        nzbb = 20;%Number of elements on the bottom side in the z direction
        nzbt = 1;%Number of elements on the top side in the z direction
        
        %computation parameters
        
        use_single = true %Determines if 32-bit (single) floats will be used instead of double size floats for the gpu calculations
        show_waitbar = true %Determines if a waitbar is shown or not
        use_gpu = true %If the gpu will be used for the matrix calculations
    end
    properties (Access=private)
        ntot = 0;%Total number of points on the mesh of the coil
        r = 0;%Will contain the a vector containing the point on the object
        dl = 0;%Will contain a vectors of the direction of the current of the points in r
        el = 0;%Will contain a normalised vecotr of the direciton of the current of the points in r
        dV = 0;%Will contain the volume related to each point in rL
        
        r_gpu = 0;%Will contain the a vector containing the points on the object, stored in gpu memory
        dl_gpu = 0;%Will contain a vectors of the direction of the current of the points in r, stored on gpu memory, only used if a coil
        el_gpu = 0;%Will contain a normalised vecotr of the direciton of the current of the points in r, stored on gpu memory, only used if a coil
        dV_gpu = 0;%Will contain the volume related to each point in rL, stored in gpu memory
    end
    methods
        function obj = set_geometry(obj,H,D1,D2,N,Rho,xpos,zpos)
            %set_geometry - Sets the parameters for the geometry of a coil
            %  set_geometry(H,D1,D2,N,Rho,xpos,zpos) makes a coil with height
            %  H,inner diameter D1, outer diameter D2, N turns, wire
            %  resistivity Rho at x position xpos and z position zpos
            obj.H = H;
            obj.D1 = D1;
            obj.D2 = D2;
            obj.N = N;
            obj.Rho = Rho;
            obj.xpos = xpos;
            obj.zpos = zpos;
        end
        
        function obj = set_nozzle_geometry(obj, H,D1,D2,Rho,xpos,zpos,taper_width,taper_angle,hole_offset,rotation)
            %set_nozzle_geometry - Sets the parameters for the geometry of a nozzle
            %  set_nozzle_geometry(H,D1,D2,N,Rho,xpos,zpos,taper_width,taper
            %  _angle,hole_offset,rotation) makes a coil with height H,inner
            %  diameter D1, outer diameter D2, at x position xpos and z
            %  position zpos, horizontal taper width taper_width, taper
            %  angle taper_angle, a hole thats offset from the center by
            %  hole_offset and i rotated by rotation degrees
            obj.H = H;
            obj.D1 = D1;
            obj.D2 = D2;
            obj.Rho = Rho;
            obj.xpos = xpos;
            obj.zpos = zpos;
            obj.taper_width = taper_width;
            obj.taper_angle = taper_angle;
            obj.hole_offset = hole_offset;
            obj.rotation = rotation;
        end
        
        function obj = set_mesh(obj, nx,nz,nphi)
            %set_mesh - Sets the parameters of a normal mesh
            %  set_mesh(nx,nz,nphi) makes a mesh with nx points in the x
            %  direction, nz point in the z direction and nphi point in the
            %  radial direction
            obj.nx = nx;
            obj.nz = nz;
            obj.nphi = nphi;
        end
        
        function obj = set_boundary_mesh(obj, nx,nz,nphi,nxBl,nxBr,nxC,nzBb,nzBt)
            %set_boundary_mesh - Sets the parameters of a normal mesh
            %  set_boundary_mesh(nx,nz,nphi,nxBl,nxBr,nxC,nzBb,nzBt) makes 
            %  a mesh with:
            %  - nx the number of elements on the bottom/top boundary mesh in the x direction
            %  - nz the number of elements on the left/right/center boundary mesh in the z direction
            %  - nphi the number of elements of the coil in radial direction
            %  - nxb the total number of elements on the left and right side in the x direction
            %  - nxbl the number of elements on the left side in the x direction
            %  - nxbr the number of elements on the right side in the x direction
            %  - nxc the number of elements on the center in the x direction
            %  - nzb the total number of elements on the top and bottom side in the z direction
            %  - nzbb the number of elements on the bottom side in the z direction
            %  - nzbt the number of elements on the top side in the z direction
            obj.nx = nx;
            obj.nz = nz;
            obj.nphi = nphi;
            obj.nxbl = nxBl;
            obj.nxbr = nxBr;
            obj.nxb = obj.nxbl+obj.nxbr;
            obj.nzbb = nzBb;
            obj.nzbt = nzBt;
            obj.nzb = obj.nzbb+obj.nzbt;
            obj.nxc = nxC;
        end
        
        

        function plot_geometry(obj,C,prefix)
            %plot_geometry - Plot the geometry of the object
            %   plot_geometry(C,prefix) plot the geometry with colour C and
            %   all direction scaled by a factor prefix
            if obj.r == 0
               error('error: cannot plot geometry. Make sure the coil has been build.') 
            end
            scatter3(obj.r(:,1)*prefix(1),obj.r(:,2)*prefix(2),obj.r(:,3)*prefix(3),C)
            daspect([1 1 1])
        end

        function plot_surface_vectors(obj,prefix)
            %plot_surface_vectors - Plots the surface vectors pointing
            %along of the current in the coil
            %   plot_surface_vectors(prefix) Plots the surface vectors 
            %   pointing in the direction of the current in the coil and
            %   all direction scaled by a factor prefix
            quiver3(obj.r(:,1)*prefix(1),obj.r(:,2)*prefix(2),obj.r(:,3)*prefix(3),obj.dl(:,1)*prefix(1),obj.dl(:,2)*prefix(2),obj.dl(:,3)*prefix(3),0 )
            daspect([1 1 1])
        end
        
        function plot_current_direction(obj)
            %plot_current_direction - Plots the normalised surface vectors
            %pointing in the direction the current in the coil
            %   plot_current_direction() Plots the normalised surface vectors
            %   pointing in the direction of the current in the coil
            quiver3(obj.r(:,1),obj.r(:,2),obj.r(:,3),obj.el(:,1),obj.el(:,2),obj.el(:,3))
            xlabel('x position (m)')
            ylabel('y position (m)')
            zlabel('z position (m)')
            daspect([1 1 1])
        end
        
        function L = calculate_self_inductance(obj)
            %calculate_self_inductance - Calculates the self inductance of
            %a coil
            %   calculate_self_inductance() Calculates the self inductance of
            %   a coil
            if obj.show_waitbar
                f = waitbar(0,'Calculating self inductance matrix');
            else
                f = 0;
            end
            O = obj.H*((obj.D2-obj.D1)/2);
            if obj.use_gpu
                M_CC = obj.inductance_matrix(obj.r,f);
            else
                M_CC = obj.inductance_matrix_gpu(obj.r_gpu,f);
            end
            L = (obj.N^2)/(obj.nx*obj.nz*O)*trace(obj.dl.'*M_CC*obj.el);
            if obj.show_waitbar
                close(f)
            end
        end
        
        function L = calculate_mutual_inductance(obj,coil)
            %calculate_mutual_inductance - Calculates the mutual inductance
            %between two coils
            %   calculate_mutual_inductance(coil) Calculates the mutual
            %   inductance between this object and a coil
            if obj.show_waitbar
                f = waitbar(0,'Calculating self inductance matrix');
            else
                f = 0;
            end
            O = obj.H*((obj.D2-obj.D1)/2);
            if obj.use_gpu
                M_CN = coil.inductance_matrix(obj.r,f);
            else
                M_CN = coil.inductance_matrix_gpu(obj.r_gpu,f);
            end
            L = (obj.N^2)/(obj.nx*obj.nz*O)*trace(obj.dl.'*M_CN*coil.el);
            if obj.show_waitbar
                close(f)
            end
        end
        
        function Z = calculate_impedance(obj,frequency)
            %calculate_impedance - Calculates the impedance of a coil
            %  calculate_impedance(frequency) Calculates the impedance of a
            %  coil at a specific frequency, without any object near it
            if obj.show_waitbar
                f = waitbar(0,'Calculating self inductance matrix');
            else
                f = 0;
            end
            O = obj.H*((obj.D2-obj.D1)/2);
            if obj.use_gpu
                M_CC = obj.inductance_matrix_gpu(obj.r_gpu,f);
            else
                M_CC = obj.inductance_matrix(obj.r,f);
            end
            Z = (obj.N^2)*1i*frequency*2*pi/(obj.nx*obj.nz*O)*trace(obj.dl.'*M_CC*obj.el)+obj.resistance();
            if obj.show_waitbar
                close(f)
            end
        end
        
        function dZ = calculate_impedance_change(obj,nozzle,frequency)
            %calculate_impedance_change - Calculates the change in
            %impedance of a coil due to an object
            %  calculate_impedance_change(nozzle, frequency) Calculates the
            %  the impedance change of a coil due to an object nozzle at a 
            %  specific frequency 
            w = frequency*2*pi;

            
            if obj.show_waitbar
                f = waitbar(0,'Calculating coupling between coil and object');
            else
                f=0;
            end
            
            if obj.use_gpu
                Mcn = nozzle.inductance_matrix_gpu(obj.r_gpu,f);
            else
                Mcn = nozzle.inductance_matrix(obj.r,f);
            end
            
            if obj.show_waitbar
                waitbar(0,f,'Calculating objects self inductance matrix');
            end
            
            if obj.use_gpu
                Mnn = nozzle.inductance_matrix_gpu(nozzle.r_gpu,f);
            else
                Mnn = nozzle.inductance_matrix(nozzle.r,f);
            end
            
            if obj.show_waitbar
                waitbar(0,f,'Calculating coupling between object and coil');
            end
            
            if obj.use_gpu
                Mnc = obj.inductance_matrix_gpu(nozzle.r_gpu,f);
            else
                Mnc = obj.inductance_matrix(nozzle.r,f);
            end
            
            if obj.show_waitbar
                waitbar(0,f,'Calculating inductance change');
            end
            
            if obj.is_circular
                O = obj.H*((obj.D2-obj.D1)/2)*pi/4;
            else
                O = obj.H*((obj.D2-obj.D1)/2);
            end
            part1 = -(obj.N^2*w^2)/(obj.nx*obj.nz*O*nozzle.Rho);
            part2 = obj.dl.'*Mcn;
            part3 = eye(nozzle.ntot)-1i*w/nozzle.Rho*Mnn;
            part4 = Mnc*obj.el;
            part5 = part3\part4;
            dZ = part1*trace(part2*part5);
            dZ = gather(dZ);
            if obj.show_waitbar
                close(f)
            end
        end

        function R = resistance(obj)
             %resistance - Calculates the DC resistance of a coil
             %   resistance() Calculates the DC resistance of this coil
             if obj.is_circular
                 A = obj.H*((obj.D2-obj.D1)/2)*pi/4/obj.N;
             else
                 A = obj.H*((obj.D2-obj.D1)/2)/obj.N;
             end
             labs = vecnorm(obj.dl,2,2);
             L = sum(labs)/obj.nx/obj.nz*obj.N;
             R = obj.Rho*L/A;
        end
        
        function depth = skin_depth(obj,f)
            %skin_depth - Analytically calculates the skin depth
            %   skin_depth(f) Calculates the skin depth with the
            %   resistivity of this object and at the frequency f
            u0 = 4*pi*1e-7;
            depth = sqrt(2*obj.Rho/u0/2/pi/f);
        end
        
        function JN_magn = calculate_induced_current_density(obj,nozzle,I,f)
            %calculate_induced_current_density - Calculates the induced
            %current density
            %   calculate_induced_current_density(nozzle,I,frequency)
            %   Calculates the current density induced in the object nozzle
            %   by this coil at current I and a specific frequency f
            w = f*2*pi;
            
            if obj.show_waitbar
                f = waitbar(0,'Calculating objects self inductance matrix');
            end
            
            if obj.use_gpu
                Mnn = nozzle.inductance_matrix_gpu(nozzle.r_gpu,f);
            else
                Mnn = nozzle.inductance_matrix(nozzle.r,f);
            end
            
            if obj.show_waitbar
                waitbar(0,f,'Calculating coupling between object and coil');
            end
            
            if obj.use_gpu
                Mnc = obj.inductance_matrix_gpu(nozzle.r_gpu,f);
            else
                Mnc = obj.inductance_matrix(nozzle.r,f);
            end  
            
            O = obj.H*((obj.D2-obj.D1)/2);
            part1 = eye(nozzle.ntot)-1i*w/nozzle.Rho*Mnn;
            part2 = 1i*w*obj.N/nozzle.Rho/O*Mnc*obj.el*I;
                        
            JN = part1\part2;
            JN = double(gather(JN));
            JN_magn = sqrt(JN(:,1).^2+JN(:,2).^2+JN(:,3).^2);
            
            if obj.show_waitbar
                close(f)
            end   
        end
        
        function plot_induced_current_density_2d(obj,nozzle,I,f,direction,prefix,offset)
            %plot_induced_current_density_2d Plot the induced current
            %density in 2d
            %   plot_induced_current_density_2d(nozzle,I,frequency,angle,
            %   prefix,offset) Plots the current density induced by this
            %   coil in the object nozzle due to a current I at a specific 
            %   frequency f in plane with direction direction. The x and z position
            %   as well as the current density all are scaled by a factor 
            %   defined by the vector prefix and offset by a factor offset.
            %   
            %   Example:
            %   plot_induced_current_density_2d(nozzle,1,1e6,0,
            %   [1e3,1e3,1e-6],[0,0,0])  Plots the induced current density in the
            %   nozzle due to a current of 1A at 1MHz running through the
            %   coil, scales the x and z axes to milli meters and the
            %   current density to MAm^-2
            JN_magn = obj.calculate_induced_current_density(nozzle,I,f);
            nozzle.plot_internal_scalar(abs(JN_magn),direction,'',prefix,offset)
        end
        
        function plot_current_density_difference(obj,nozzle,I,f,Y,X,J,direction,prefix,offset)
            %plot_current_density_difference - Plot the difference with an
            %externally calculated current density
            %   plot_current_density_difference(nozzle,I,frequency,Y,X,J,
            %   direction,prefix,offset) Plots the difference between the 
            %   current density induced by this coil in the object nozzle
            %   due to a current I at a specific frequency f in plane with 
            %   angle direction and an externally calculated current
            %   density J at position Y and X. The x and z position as well
            %   as the current density all are scaled by a factor defined 
            %   by the vector prefix and offset by a factor offset.
            
            n = round(((unwrap(direction)+pi)/2/pi)*(nozzle.nphi-1))+1;
            use = n:nozzle.nphi:nozzle.ntot;
            magn = sqrt(nozzle.r(use,1).^2+nozzle.r(use,2).^2);
            z = nozzle.r(use,3);
            
            JN_magn = obj.calculate_induced_current_density(nozzle,I,f);
            
            J_intp = interp2(Y, X, abs(J),z,magn);
            
            use2 = ~isnan(J_intp);
            
            J_diff = abs(JN_magn(use))-J_intp;
            
            obj.plot_trisurf(prefix(1)*magn(use2)+offset(1),prefix(2)*z(use2)+offset(2),prefix(3)*J_diff(use2)+offset(3));
            
        end
        
        function plot_current_density_difference_ratio(obj,nozzle,I,f,Y,X,J,direction,prefix,offset)
            %plot_current_density_difference_ratio - Plot the difference with an
            %externally calculated current density
            %   plot_current_density_difference_ratio(nozzle,I,frequency,Y,X,J,
            %   direction,prefix,offset) Plots the relative error between the 
            %   current density induced by this coil in the object nozzle
            %   due to a current I at a specific frequency f in plane with 
            %   angle direction and an externally calculated current
            %   density J at position Y and X. The x and z position as well
            %   as the current density all are scaled by a factor defined 
            %   by the vector prefix and offset by a factor offset. The
            %   error is relative to the externally calculated current
            %   density
            n = round(((unwrap(direction)+pi)/2/pi)*(nozzle.nphi-1))+1;
            use = n:nozzle.nphi:nozzle.ntot;
            magn = sqrt(nozzle.r(use,1).^2+nozzle.r(use,2).^2);
            z = nozzle.r(use,3);
            
            JN_magn = obj.calculate_induced_current_density(nozzle,I,f);
            
            J_intp = interp2(Y, X, abs(J),z,magn);
            
            use2 = ~isnan(J_intp);
            
            J_diff = abs(abs(JN_magn(use))-J_intp)./J_intp*100;
            
            obj.plot_trisurf(prefix(1)*magn(use2)+offset(1),prefix(2)*z(use2)+offset(2),prefix(3)*J_diff(use2)+offset(3));
        end
        
        function plot_trisurf(~,x,y,z)
            %plot_trisurf - Plot a 1D vector of values in 2D
            %   plot_trisurf(x,y,z) Makes a 2D plot in which z determines
            %   the colour and x and y determine the position of that
            %   colour
            dt = delaunayTriangulation(x,y) ;
            tri = dt.ConnectivityList ;
            xi = dt.Points(:,1) ; 
            yi = dt.Points(:,2) ; 
            F = scatteredInterpolant(x,y,z);
            zi = F(xi,yi) ;
            trisurf(tri,xi,yi,zi) 
            view(2)
            shading interp
        end 
        
            
        function plot_internal_scalar(obj,phi,direction,label,prefix,offset)
            %plot_internal_scalar - Plot a scalar from inside this object
            %in 2D
            %   plot_internal_scalar(phi,direction,label,prefix,offset)
            %   Plots the variable Phi in the plane with direction direction.
            %   It gives the c axis the label label and scales x,z,Phi with 
            %   the values of the vector prefix and gives them an offset from
            %   the vector offset.
            n = round(((unwrap(direction)+pi)/2/pi)*(obj.nphi-1))+1;
            use = n:obj.nphi:obj.ntot;

            magn = sqrt(obj.r(use,1).^2+obj.r(use,2).^2);
            z = obj.r(use,3);
            obj.plot_trisurf(magn*prefix(1)+offset(1),z*prefix(2)+offset(2),real(phi(use))*prefix(3)+offset(3));
            hold on
            scatter3(magn*prefix(1)+offset(1), z*prefix(2)+offset(2),max(phi(use))*ones(length(use),1)*prefix(3)+offset(3),8,'r','filled')
            c = colorbar;
            c.Label.String = label;
        end
        
        function L = inductance_matrix(obj,rd,f)
            %inductance_matrix - Calculate the summation of currents matrix
            %   inductance_matrix(rd,f) Calculates the summation of
            %   currents matrix to relate current densities in the current
            %   object to the A field at the locations specified in rd
            u0 = 4*pi*1e-7;
            nd = length(rd);
            L = zeros(nd,obj.ntot);
            parfor i1 = 1:nd
                rd_m = ones(obj.ntot,1)*rd(i1,:);
                r_dif = rd_m-obj.r;
                r_abs = vecnorm(r_dif,2,2);   
                L(i1,:) = u0/4/pi./r_abs.*obj.dV;
            end
            L(isinf(L)) = 0;
        end
        
        function L = inductance_matrix_gpu(obj,rd,f)
            %inductance_matrix_gpu - Calculate the summation of currents
            %matrix using the gpu
            %   inductance_matrix(rd,f) Calculates the summation of
            %   currents matrix to relate current densities in the current
            %   object to the A field at the locations specified in rd inside
            %   using gpu
            if obj.use_single
                u0 = gpuArray(single(4*pi*1e-7));
                nd = length(rd);
                L = gpuArray(single(zeros(nd,obj.ntot)));
                for i1 = 1:nd
                    if obj.show_waitbar
                        if mod(i1,100) == 0
                            waitbar(i1/nd,f)
                        end
                    end
                    rd_m = gpuArray(single(ones(obj.ntot,1)))*rd(i1,:);
                    r_dif = rd_m-obj.r_gpu;
                    r_abs = vecnorm(r_dif,2,2);   
                    L(i1,:) = u0/4/pi./r_abs.*obj.dV_gpu;
                end
            else
                u0 = gpuArray(4*pi*1e-7);
                nd = length(rd);
                L = gpuArray(zeros(nd,obj.ntot));
                for i1 = 1:nd
                    if obj.show_waitbar
                         if mod(i1,100) == 0
                            waitbar(i1/nd,f)
                         end
                    end
                    rd_m = gpuArray(ones(obj.ntot,1))*rd(i1,:);
                    r_dif = rd_m-obj.r_gpu;
                    r_abs = vecnorm(r_dif,2,2);   
                    L(i1,:) = u0/4/pi./r_abs.*obj.dV_gpu;
                end
            end
            L(isinf(L)) = 0;
        end
        
        function obj = build_coil(obj)
            %build_coil - Builds the mesh of a coil
            %   build_coil() Builds the mesh of a coil. 
            if obj.show_waitbar
                f = waitbar(0,'Building coil');
            end
            
            obj.is_coil = 1;
            obj.is_circular = 0;
            
            obj.ntot = obj.nz*obj.nx*obj.nphi;
            obj.r = zeros(obj.ntot,3);
            obj.dl = zeros(obj.ntot,3);
            obj.dV = zeros(obj.ntot,1);

            dx = (obj.D2-obj.D1)/2/obj.nx;
            x = obj.D1/2+dx/2:dx:obj.D2/2;

            dphi = 2*pi/obj.nphi;
            phi = dphi/2:dphi:2*pi;

            dz = obj.H/obj.nz;
            z = dz/2:dz:obj.H-dz/2;

            for i1 = 1:obj.nz
                if obj.show_waitbar
                    waitbar(i1/obj.nz,f)
                end
                for i2 = 1:obj.nx
                    for i3 = 1:obj.nphi
                        loc = (i1-1)*obj.nx*obj.nphi+(i2-1)*obj.nphi+i3;
                        obj.r(loc,1) = x(i2)*cos(phi(i3))+obj.xpos;
                        obj.r(loc,2) = x(i2)*sin(phi(i3));
                        obj.r(loc,3) = z(i1)+obj.zpos;
                        
                        obj.dl(loc,1) = -x(i2)*sin(phi(i3))*dphi;
                        obj.dl(loc,2) = x(i2)*cos(phi(i3))*dphi;
                        obj.dl(loc,3) = 0;
                        
                        obj.el(loc,1) = -sin(phi(i3));
                        obj.el(loc,2) = cos(phi(i3));
                        obj.el(loc,3) = 0;
                        
                        dS = ((x(i2)+dx/2)^2-(x(i2)-dx/2)^2)*pi/obj.nphi;
                        obj.dV(loc) = dz*dS;
                    end
                end
            end
            
            if obj.use_single 
                obj.r_gpu = gpuArray(single(obj.r));
                obj.dl_gpu = gpuArray(single(obj.dl));
                obj.el_gpu = gpuArray(single(obj.el));
                obj.dV_gpu = gpuArray(single(obj.dV));
            else
                obj.r_gpu = gpuArray(obj.r);
                obj.dl_gpu = gpuArray(obj.dl);
                obj.el_gpu = gpuArray(obj.el);
                obj.dV_gpu = gpuArray(obj.dV);
            end
            if obj.show_waitbar
                close(f)
            end
        end
        
        function obj = build_single_loop(obj)
            %build_single_loop - Builds the mesh of a single loop
            %   build_single_loop() Builds the mesh of a single loop, which is
            %   different from a coil in that the wire is approximately 
            %   circular instead of rectangular
            if obj.show_waitbar
                f = waitbar(0,'Building coil');
            end
            
            obj.is_coil = 1;
            obj.is_circular = 1;
            
            %calculate coil position vector
            obj.ntot = obj.nz*obj.nx*obj.nphi;
            obj.r = zeros(obj.ntot,3);
            obj.dl = zeros(obj.ntot,3);
            obj.dV = zeros(obj.ntot,1);

            r_circ = (obj.D2-obj.D1)/4;
            r_avg = (obj.D2+obj.D1)/4;
            

            dphi = 2*pi/obj.nphi;
            phi = dphi/2:dphi:2*pi;

            dz = obj.H/obj.nz;
            z = -obj.H/2+dz/2:dz:obj.H/2-dz/2;

            for i1 = 1:obj.nz
                if obj.show_waitbar
                    waitbar(i1/obj.nz,f)
                end
                y_circ = z(i1);
                x_circ = sqrt(r_circ^2-y_circ^2);
                ddxmax = y_circ/x_circ*dz;
                
                dx = 2*x_circ/obj.nx;
                x = r_avg-x_circ+dx/2:dx:r_avg+x_circ-dx/2;
                for i2 = 1:obj.nx
                    ddx1 = ddxmax*x(i2)-dx/2/x_circ;
                    ddx2 = ddxmax*x(i2)+dx/2/x_circ;
                    for i3 = 1:obj.nphi
                        loc = (i1-1)*obj.nx*obj.nphi+(i2-1)*obj.nphi+i3;
                        obj.r(loc,1) = x(i2)*cos(phi(i3))+obj.xpos;
                        obj.r(loc,2) = x(i2)*sin(phi(i3));
                        obj.r(loc,3) = z(i1)+obj.zpos;
                        
                        obj.dl(loc,1) = -x(i2)*sin(phi(i3))*dphi;
                        obj.dl(loc,2) = x(i2)*cos(phi(i3))*dphi;
                        obj.dl(loc,3) = 0;
                        
                        obj.el(loc,1) = -sin(phi(i3));
                        obj.el(loc,2) = cos(phi(i3));
                        obj.el(loc,3) = 0;
                        
                        xl1 = x(i2)-dx/2-ddx1;
                        xl2 = x(i2)-dx/2+ddx1;
                        xr1 = x(i2)+dx/2-ddx2;
                        xr2 = x(i2)+dx/2+ddx2;
                        obj.dV(loc) = 1/3*pi*(xr1^2+xr1*xr2+xr2^2-xl1^2-xl1*xl2-xl2^2)*dz/obj.nphi;
                    end
                end
            end
            
            if obj.use_single 
                obj.r_gpu = gpuArray(single(obj.r));
                obj.dl_gpu = gpuArray(single(obj.dl));
                obj.el_gpu = gpuArray(single(obj.el));
                obj.dV_gpu = gpuArray(single(obj.dV));
            else
                obj.r_gpu = gpuArray(obj.r);
                obj.dl_gpu = gpuArray(obj.dl);
                obj.el_gpu = gpuArray(obj.el);
                obj.dV_gpu = gpuArray(obj.dV);
            end
            
            if obj.show_waitbar
                close(f)
            end
        end
        
        function obj = build_nozzle_surface(obj,depth)
            %build_nozzle_surface - Builds the mesh of a nozzle
            %   build_nozzle_surface(depth) Builds the boundary mesh of a 
            %   nozzle with the boundary mesh extending up to a depth depth
            if obj.show_waitbar
                f = waitbar(0,'Building surface coil');
            end
            
            obj.is_coil = 0;
            obj.is_circular = 0;
            %calculate coil position vector
            obj.ntot = obj.nzb*obj.nx*obj.nphi+obj.nz*obj.nxb*obj.nphi+obj.nz*obj.nxc*obj.nphi;
            obj.r = zeros(obj.ntot,3);
            obj.dl = zeros(obj.ntot,3);
            obj.dV = zeros(obj.ntot,1);
            
            dphi = 2*pi/obj.nphi;
            phi = dphi/2:dphi:2*pi;
            
            dzb = 2*depth/obj.nzb;
            n1 = obj.nzbb;
            n2 = obj.nzbt;
            z1 = dzb/2:dzb:dzb*(n1-0.5);
            z2 = obj.H-dzb*(n2-0.5):dzb:obj.H;
            zb = [z1,z2];
            
            dz = (obj.H-2*depth)/obj.nz;
            z = dzb*n1+dz/2:dz:obj.H-dzb*n2-dz/2;
               
            for i3 = 1:obj.nphi
                if phi(i3) <= pi  
                    xstart(i3) = obj.hole_offset*sin(phi(i3))+sqrt((obj.D1/2)^2 - (obj.hole_offset*cos(phi(i3)))^2);
                else
                    xstart(i3) = obj.hole_offset*sin(phi(i3))+sqrt((obj.D1/2)^2 - (obj.hole_offset*cos(2*pi-phi(i3)))^2);
                end
            end
            
            %calculate x locations for the mesh on the left and right side of the
            %nozzle mesh
            for i2 = 1:obj.nz
                if obj.show_waitbar
                    waitbar(i2/obj.nz/5,f)
                end
                for i3 = 1:obj.nphi
                    taper_height = obj.taper_width*tan(obj.taper_angle);
                    
                    dxb = 2*depth/obj.nxb;
                    
                    n1 = obj.nxbl;
                    n2 = obj.nxbr;
                    
                    if z(i2) > taper_height
                        dxc{i2,i3} = (obj.D2/2-xstart(i3)-2*depth)/obj.nxc;
                        if dxc{i2,i3} < 0
                            error("warning: can't build boundary mesh, skind depth too large")
                        end
                        
                        x1 = xstart(i3)+dxb/2:dxb:xstart(i3)+dxb*(n1-0.5);
                        x2 = obj.D2/2-dxb*(n2-0.5):dxb:obj.D2/2-dxb/2;
                        xb{i2,i3} = [x1,x2];
                        xc{i2,i3} = xstart(i3)+dxb*n1+dxc{i2,i3}/2:dxc{i2,i3}:obj.D2/2-dxb*(n2-0.5)-dxc{i2,i3}/2;
                    else
                        xstop = obj.D2/2-obj.taper_width/taper_height*(taper_height-z(i2));
                        
                        dxc{i2,i3} = (xstop-xstart(i3)-2*depth)/obj.nxc;
                        
                        x1 = xstart(i3)+dxb/2:dxb:xstart(i3)+dxb*(n1-0.5);
                        x2 = xstop-dxb*(n2-0.5):dxb:xstop;
                        xb{i2,i3} = [x1,x2];
                        xc{i2,i3} = xstart(i3)+dxb*n1+dxc{i2,i3}/2:dxc{i2,i3}:xstop-dxb*n2-dxc{i2,i3}/2;
                    end
                end
            end
            
            %calculate r,dl,el,ds,dV for the left and right side of the
            %mesh
            for i1 = 1:obj.nxb
                if obj.show_waitbar
                    waitbar(i1/obj.nxb/5+0.2,f)
                end
                for i2 = 1:obj.nz
                    for i3 = 1:obj.nphi
                        loc = (i1-1)*obj.nz*obj.nphi+(i2-1)*obj.nphi+i3;
                        obj.r(loc,1) = xb{i2,i3}(i1)*cos(phi(i3));
                        obj.r(loc,2) = xb{i2,i3}(i1)*sin(phi(i3));
                        obj.r(loc,3) = z(i2);
                        
                        if z(i2) >= taper_height
                            ddx = 0;
                        else
                            ddx = obj.taper_width/taper_height*dz;
                        end
                        xl1 = xb{i2,i3}(i1)-dxb/2-ddx;
                        xl2 = xb{i2,i3}(i1)-dxb/2+ddx;
                        xr1 = xb{i2,i3}(i1)+dxb/2-ddx;
                        xr2 = xb{i2,i3}(i1)+dxb/2+ddx;
                        obj.dV(loc) = 1/3*pi*(xr1^2+xr1*xr2+xr2^2-xl1^2-xl1*xl2-xl2^2)*dz/obj.nphi;
                    end
                end
            end
            

            for i1 = 1:obj.nxc
                if obj.show_waitbar
                    waitbar(i1/obj.nxc/5+0.4,f)
                end
                for i2 = 1:obj.nz
                    for i3 = 1:obj.nphi
                        loc = (i1-1)*obj.nz*obj.nphi+(i2-1)*obj.nphi+i3+obj.nz*obj.nxb*obj.nphi;
                        obj.r(loc,1) = xc{i2,i3}(i1)*cos(phi(i3));
                        obj.r(loc,2) = xc{i2,i3}(i1)*sin(phi(i3));
                        obj.r(loc,3) = z(i2);
                        
                        if z(i2) >= taper_height
                            ddx = 0;
                        else
                            xstop = obj.D2/2-obj.taper_width/taper_height*(taper_height-z(i2));
                            xcstop = xstop-dxb*n2;
                            xcstart = xstart(i3)+dxb*n1;
                            ddx = obj.taper_width/taper_height*dz*(xc{i2,i3}(i1)-xcstart)/(xcstop-xcstart);
                        end
                        
                        xl1 = xc{i2,i3}(i1)-dxc{i2,i3}/2-ddx;
                        xl2 = xc{i2,i3}(i1)-dxc{i2,i3}/2+ddx;
                        xr1 = xc{i2,i3}(i1)+dxc{i2,i3}/2-ddx;
                        xr2 = xc{i2,i3}(i1)+dxc{i2,i3}/2+ddx;
                        obj.dV(loc) = 1/3*pi*(xr1^2+xr1*xr2+xr2^2-xl1^2-xl1*xl2-xl2^2)*dz/obj.nphi;
                    end
                end
            end
            
            
            
            %calculate x locations for the mesh on the top and bottom side of the
            %nozzle mesh
            for i3 = 1:obj.nphi
                if phi(i3) <= pi  
                    xstart(i3) = obj.hole_offset*sin(phi(i3))+sqrt((obj.D1/2)^2 - (obj.hole_offset*cos(phi(i3)))^2);
                else
                    xstart(i3) = obj.hole_offset*sin(phi(i3))+sqrt((obj.D1/2)^2 - (obj.hole_offset*cos(2*pi-phi(i3)))^2);
                end
            end
            
            for i1 = 1:obj.nzb
                if obj.show_waitbar
                    waitbar(i1/obj.nzb/5+0.5,f)
                end
                for i3 = 1:obj.nphi
                    taper_height = obj.taper_width*tan(obj.taper_angle);
 
                    if zb(i1) >= taper_height

                        dx{i1,i3} = (obj.D2/2-xstart(i3))/obj.nx;
                        x{i1,i3} = xstart(i3)+dx{i1,i3}/2:dx{i1,i3}:obj.D2/2-dx{i1,i3}/2;
                        
                    else
                        xstop = obj.D2/2-obj.taper_width/taper_height*(taper_height-zb(i1));

                        dx{i1,i3} = (xstop-xstart(i3))/obj.nx;
                        x{i1,i3} = xstart(i3)+dx{i1,i3}/2:dx{i1,i3}:xstop;
                        
                        
                    end
                end
            end
            
            %calculate r,dl,el,ds,dV for the top and bottom side of the
            %mesh
            for i1 = 1:obj.nzb
                if obj.show_waitbar
                    waitbar(i1/obj.nzb/5+0.8,f)
                end
                for i2 = 1:obj.nx
                    for i3 = 1:obj.nphi
                        loc = (i1-1)*obj.nx*obj.nphi+(i2-1)*obj.nphi+i3+obj.nz*obj.nxb*obj.nphi+obj.nz*obj.nxc*obj.nphi;
                        obj.r(loc,1) = x{i1,i3}(i2)*cos(phi(i3));
                        obj.r(loc,2) = x{i1,i3}(i2)*sin(phi(i3));
                        obj.r(loc,3) = zb(i1);
                        
                        if zb(i1) >= taper_height
                            ddx = 0;
                        else
                            xstop = obj.D2/2-obj.taper_width/taper_height*(taper_height-zb(i1));
                            ddx = obj.taper_width/taper_height*dzb*(x{i1,i3}(i2)-xstart(i3))/(xstop-xstart(i3));
                        end
                        xl1 = x{i1,i3}(i2)-dx{i1,i3}/2-ddx;
                        xl2 = x{i1,i3}(i2)-dx{i1,i3}/2+ddx;
                        xr1 = x{i1,i3}(i2)+dx{i1,i3}/2-ddx;
                        xr2 = x{i1,i3}(i2)+dx{i1,i3}/2+ddx;
                        z1 = zb(i1)-dzb/2;
                        z2 = zb(i1)+dzb/2;
                        obj.dV(loc) = 1/3*pi*(xr1^2+xr1*xr2+xr2^2-xl1^2-xl1*xl2-xl2^2)*dzb/obj.nphi;
                        
                    end
                end
            end
            

            rtemp = obj.r(:,1);
            obj.r(:,1) = obj.r(:,2);
            obj.r(:,2) = rtemp;
            
            rotation_rad =  obj.rotation*pi/180;
            obj.r(:,1) = cos(rotation_rad)*obj.r(:,1)-sin(rotation_rad)*obj.r(:,3);
            obj.r(:,3) = sin(rotation_rad)*obj.r(:,1)+cos(rotation_rad)*obj.r(:,3);
            
            obj.r(:,1) = obj.r(:,1) + obj.xpos;
            obj.r(:,3) = obj.r(:,3) + obj.zpos;
            
            if obj.use_single 
                obj.r_gpu = gpuArray(single(obj.r));
                obj.dV_gpu = gpuArray(single(obj.dV));
            else
                obj.r_gpu = gpuArray(obj.r);
                obj.dV_gpu = gpuArray(obj.dV);
            end
            if obj.show_waitbar
                close(f)
            end
        end
    end
end