import numpy as np

def adjust_vz(vz, x_coords, y_coords, ref_x, ref_y): #set reference point to zero, make potential positive l.281

    ix = np.argmin(np.abs(x_coords - ref_x))
    iy = np.argmin(np.abs(y_coords - ref_y))

    reference_value = vz[ix, iy]

    return -(vz - reference_value)
    

def build_supercell(vz, cells_x, cells_y): # here the adjusted vz , l.308 onwards
    periodic_part = vz[:-1, :-1]

    repeated = np.tile(periodic_part, (cells_x, cells_y))

    # Restore one periodic boundary point on each side.
    vz_sc = np.empty(
        (repeated.shape[0] + 1, repeated.shape[1] + 1),
        dtype=repeated.dtype,
    )

    vz_sc[:-1, :-1] = repeated
    vz_sc[-1, :-1] = repeated[0, :]
    vz_sc[:-1, -1] = repeated[:, 0]
    vz_sc[-1, -1] = repeated[0, 0]

    return vz_sc

def make_real_space_grids(nx, ny, dx, dy, nx_sc, ny_sc):
    xsc_npt = nx_sc*(nx-1)+1 # l.324
    ysc_npt = ny_sc*(ny-1)+1
    x = np.arange(xsc_npt, dtype=float) * dx #l. 340
    y = np.arange(ysc_npt, dtype=float) * dy 
    return x, y

def make_reciprocal_space_grid():
    # l.399 
    
      do i=1,ngx
     kxgrid(i)=-real((ngx/2),kind=dp)*dkx+real((i-1),kind=dp)*dkx
  end do
  do j=1,ngy
     kygrid(j)=-real((ngy/2),kind=dp)*dky+real((j-1),kind=dp)*dky
  end do
  do i=1,ngx
     do j=1,ngy
        kr(i,j)=sqrt(kxgrid(i)**2 + kygrid(j)**2)
     end do
  end do
    pass

def fft_forward(vz, dx, dy):
    return (
        dx
        * dy
        * np.fft.fftshift(np.fft.fft2(vz))
    )


def fft_backward(vz_k, dx, dy):
    return (
        np.fft.ifft2(np.fft.ifftshift(vz_k))
        / (dx * dy)
    )






