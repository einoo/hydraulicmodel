# hydraulicmodel
Flood modeling by 2D shallow water equation. 
Refer to Hunter et al (2005), Bates et al. (2010). 

Please refer to the following paper on the detailed implementation and results. 

```
@Article{ouyang2022employment,
  author  = {Mao Ouyang and Shunji Kotsuki and Yuka Ito and Tomochika Tokunaga},
  journal = {Journal of Hydrology: Regional Studies},
  title   = {Employment of hydraulic model and social media data for flood hazard assessment in an urban city},
  year    = {2022},
  issn    = {2214-5818},
  pages   = {101261},
  volume  = {44},
  doi     = {10.1016/j.ejrh.2022.101261},
}
```

## Diffusive wave approximation

<img src="https://render.githubusercontent.com/render/math?math=\frac{\partial h ^{i, j}}{\partial t} = \frac{Q_x^{i-1, j} - Q_x^{i, j} + Q_y^{i, j-1} - Q_y^{i, j}}{\Delta x \Delta y}">

<img src="https://render.githubusercontent.com/render/math?math=Q_x^{i, j} = \frac{h_{flow}^{5/3}}{n} \left( \frac{h^{i-1, j} - h^{i, j}}{\Delta x} \right) ^{1/2} \Delta y">

<img src="https://render.githubusercontent.com/render/math?math=Q_y^{i, j} = \frac{h_{flow}^{5/3}}{n} \left( \frac{h^{i-1, j} - h^{i, j}}{\Delta y} \right) ^{1/2} \Delta x">

## Local inertail approaximation

<img src="https://render.githubusercontent.com/render/math?math=\frac{\partial h ^{i, j}}{\partial t} = \frac{Q_x^{i-1, j} - Q_x^{i, j} + Q_y^{i, j-1} - Q_y^{i, j}}{\Delta x \Delta y}">

<img src="https://render.githubusercontent.com/render/math?math=q_{t + \Delta t} = \frac{q_t - gh_{flow, t} \Delta t \frac{\partial h_t}{\partial x}} {1 + g n^2 \Delta t q_t / h_{flow, t}^{7/3}}">


## Time steps in two methods

- Diffusive wave approximation
<img src="https://render.githubusercontent.com/render/math?math=\Delta t = \frac{\Delta x^2}{4} min \left( \frac{2n}{h_{flow}^{5/3}}|\frac{\partial h}{\partial x}|^{1/2}, \frac{2n}{h_{flow}^{5/3}}|\frac{\partial h}{\partial y}|^{1/2} \right)">

- Local inertial approximation
<img src="https://render.githubusercontent.com/render/math?math=\Delta t = \alpha \frac{\Delta x}{\sqrt{g (h-z)_{max}}}">

## Implementation of these two approximations for simulating the flood in Mobara city on October 25, 2019

Setting of the domain and boundary conditions
![domain](./fig/domain.png)

Results of the maximum inundation area calculated by the diffusive wave approximation, compared with the flood extent reported by Chiba Office
![diffusive wave](./fig/dw.png)

Results of the inundation range considering the flood recession by the local inertial approximation, compared with the flood extent estimated by Geospatial Information Authority of Japan (GSI). 
![local inertial](./fig/li.png)

## Evolution of flood volume
![flood volume](./fig/fv.png)
