# High synaptic threshold for NMDA spike generation in human layer 2/3 pyramidal neurons

To run:
	1. Compile mod files in the folder into a new "special". Alternatively, you can drag and drop the folder of interest into the mknrndll application icon in NEURON.
	2. Run Terminal and change to directory containing this file.
	3. Execute the following:
	mod/x86_64/special init-human-simulation.hoc, init-mouse-simulation.hoc, 	ball_and_stick_model.hoc.
	4. In Figure 3&4 folder, you’ll find two folders: Human and Mouse simulations. Compile the mod 	files 	into the respective folders and run init-human/mouse-simulation.hoc file. Change the 	number of synapses from the parameter section in the NEURON GUI. To overlay the plots 	for different synapse number-left click on the simulation window of soma voltage panel 	select ‘view plot’ and 	then turn on ’keep lines’. 
	5. Choose different basal dendrites from the pploc section list. For human simulations: {dend [55], dend [29], dend [18], dend [26], dend [24], dend [46]}.  
	For mouse simulations: { a1_1111, a1_1211, a1_1121, a4_1111, a2_1112}.
	6.For Synaptic location simulations in Figure 4, set number of synapses, N to 1 and then vary gmax between (5 to 35 nS). The synaptic location value is a normalized position for the synapse point process on a 	specific dendritic branch chosen from the section list.
	7. In Figure 5 folder, you will be able to run ball_and_stick_model.hoc file after mod files compilation.
	vary the parameter of your choice from the NEURON GUI and repeat step 4 to view overlaid traces.
