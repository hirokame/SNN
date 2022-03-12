# Spiking Neural Network

Simulation of Basal ganglia circuit.  

Network
- Contains Cortex, Striatum, GPe, SNr, Thalamus.  
- Cortex have an excitatory projections to the striatum which has the two sutypes of neuron(D1 and D2 expressing neuron).  
- Striatum D1 expressing neuron have an direct inhibitory projections to SNr(direct pathway).  
- Straitum D2 expressing neuron have an indirect inhibitory projections to SNr(indirect pathway) via GPe region.  
- SNr projects to Thalamus inhibitory, then Thalamus projects back to Cortex excitatory.  

Neuron
- Each neurons are implemented in integrate and Fire neuron.  
- time step = 0.1ms
- refractory period = 2ms
- resting potencial = -70mV
- threshold = -50mV
