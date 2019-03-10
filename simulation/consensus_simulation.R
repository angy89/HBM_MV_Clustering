# The best way to run this code is to let R execute it as a whole by typing
# source("ConsensusSimulation.R") 
# at the R command prompt. The script will run the simulation code, producing a lot of diagnostic output
# and some figures demoinstrating what it is doing. At the end of the simulation it will save the
# simulation results and stop at the 
# stop("Simulation is now done.");
# command (around line 510). The rest (statistical analysis of the results) is best run by copy-pasting
# the code into the R session.

#setwd("C:/Documents and Settings/plangfelder/My Documents/Work/GeneExpressionSimulation/");
source("NetworkFunctions-Simulation.R");

# Set overall parameters
options(stringsAsFactors = FALSE);
set.seed(12);
No.Samples = 120;
Set.Labels = c("Set 1", "Set 2");
Subsets = c(rep(1, times = 55), rep(2, times = 65));
No.Sets = length(Set.Labels);

No.Genes = 3000;
MinSimModSize = 50;
MinProportion = MinSimModSize/No.Genes;
GreyProportion = 0.05;

No.SimulationLevels = 6;

BackgrNoise = seq(from = 0, to = 0.4, length.out = No.SimulationLevels);
NOrderedSmallLayers = as.integer(seq(from = 0.5, to = 4.5, length.out = No.SimulationLevels));
NUnorderedSmallLayers = as.integer(seq(from = 0.5, to = 6.2, length.out = No.SimulationLevels));
AverageNGenesInSmallModule = as.integer(seq(from = 10, to = 25, length.out = No.SimulationLevels));
AverageExprInSmallModule = seq(from = 0.2, to = 0.60, length.out = No.SimulationLevels);

No.Simulations = 100;

DetBranchCuts = seq(from = 0.965, to = 0.995, length.out = 4);
No.DetBranchCuts = length(DetBranchCuts);
SetBranchCuts = seq(from = 0.96, to = 0.99, length.out = 3);
No.SetBranchCuts = length(SetBranchCuts);

MaxModules = 21;

PlotDir = "";
PlotFileBase = "Simulation-";
PlotsPerLevel = 3;

# Initialize simulation statistics...

SimStats = list(SimModIsConsensus = NULL, SimModPValues = NULL, SimModNextPValues = NULL,
                HighestPValIsGrey = NULL,
                SimModPropInClosest = NULL, SimModPropInGrey = NULL, SimModPropNotInGreyOrClosest = NULL, 
                Seed_ClosestPC_Cor = NULL, PCCor = NULL, 
                ConnMax = NULL, ConnAssdAve = NULL, ConnUnassdAve = NULL);

#SimModIsConsensus : Simulated module is consensus/simulated in set 1/simulated in set 2 : TRUE or FALSE
#SimModPValues : Simulated module's highest p-value in detected modules (inc. grey)
#SimModNextPValues: Simulated module's next-highest p-value in detected modules (incl. grey)
#HighestPValIsGrey: True or false
#SetHighestPValIsGrey: T or F: for each set separately perform module detection and 
#      simulated <--> detected cross assignment.
# ConnMax: maximum connectivity of genes in each _simulated_ module
# ConnAssdAve: average connectiviy of all genes in a given simulated module that are assigned to a
#          detected consensus module
# ConnUnassdAve: average connectiviy of all genes in a given simulated module that are detected consensus
#          grey 

SimStats$SimModIsConsensus = array(NA,  dim = c(MaxModules, No.Sets+1, No.Simulations, No.SimulationLevels));
SimStats$SimModPValues = array(NA, dim = c(MaxModules, No.Sets+1, No.DetBranchCuts, No.Simulations, 
                                           No.SimulationLevels));
SimStats$SimModNextPValues = array(NA, dim = c(MaxModules, No.Sets+1, No.DetBranchCuts, No.Simulations, 
                                               No.SimulationLevels));
SimStats$HighestPValIsGrey = array(NA, dim = c(MaxModules, No.Sets+1, No.DetBranchCuts, No.Simulations, 
                                               No.SimulationLevels));
SimStats$SimModPropInClosest = array(NA, dim = c(MaxModules, No.Sets+1, No.DetBranchCuts, No.Simulations, 
                                                 No.SimulationLevels));
SimStats$SimModPropInGrey = array(NA, dim =c(MaxModules, No.Sets+1, No.DetBranchCuts, No.Simulations, 
                                             No.SimulationLevels));
SimStats$SimModPropNotInGreyOrClosest = array(NA, dim =c(MaxModules, No.Sets+1, No.DetBranchCuts, 
                                                         No.Simulations, No.SimulationLevels));
SimStats$Seed_ClosestPC_Cor = array(NA, dim =c(MaxModules, No.Sets+1, No.DetBranchCuts, No.Simulations, 
                                               No.SimulationLevels));
SimStats$PCCor = vector(mode = "list", 
                        length = No.Sets * No.DetBranchCuts * No.Simulations * No.SimulationLevels);
dim(SimStats$PCCor) = c(No.Sets, No.DetBranchCuts, No.Simulations, No.SimulationLevels);

SimStats$ConnMax = array(NA, dim = c(MaxModules, No.Sets, No.DetBranchCuts, No.Simulations, 
                                     No.SimulationLevels));
SimStats$ConnAssdAve = array(NA, dim = c(MaxModules, No.Sets, No.DetBranchCuts, No.Simulations, 
                                         No.SimulationLevels));
SimStats$ConnUnassdAve = array(NA, dim = c(MaxModules, No.Sets, No.DetBranchCuts, No.Simulations, 
                                           No.SimulationLevels));

SimulationLevel = 1; simulation = 1;

for (SimulationLevel in 1:No.SimulationLevels) 
{
  
  for (simulation in 1:No.Simulations)
  {
    set.seed(No.Simulations*SimulationLevel + simulation);
    # Create the causal network: anchors and causal matrix
    
    NModuleBlocks = 5;
    MinNModsInBlock = 2;
    MaxNModsInBlock = 5;
    
    # Generate number of modules in each block, but no more than MaxModules
    
    No.Nodes = MaxModules + 2;
    print.flush("Generating module block membership...");
    while (No.Nodes>=MaxModules)  # MaxModules must also accomodate grey, which is not counted in here.
    {
      NModsInBlock = as.integer(runif(NModuleBlocks, min = MinNModsInBlock, max = MaxNModsInBlock));
      No.Nodes = sum(NModsInBlock);
    }
    
    # Generate one anchor for each block
    
    AnchorInd = c(1, cumsum(NModsInBlock[1:(NModuleBlocks-1)])+1);
    
    No.Anchors = length(AnchorInd);
    
    AnchorVectors = matrix( rnorm(n=No.Samples*No.Anchors), nrow = No.Samples, ncol = No.Anchors);
    
    # Generate the causal dependence matrix
    
    Noise = 1; MeanCause = 0.2; SigmaCause = 0.5;
    MinTotCause = 0.1;
    Cause = matrix(0, nrow = No.Nodes, ncol = No.Nodes);
    print.flush("Generating module block causal matrices...");
    for (block in 1:NModuleBlocks)
    {
      BlockCause = matrix(0, nrow = NModsInBlock[block], ncol = NModsInBlock[block]);
      nentries = sum(lower.tri(BlockCause))
      # This condition makes sure that within each block, the sum of causes for each vector is at least MTC
      while (min(apply(abs(BlockCause), 1, sum)) < MinTotCause)
      {
        # To satisfy the minimum total cause condition, set the last diagonal element
        BlockCause[1,1] = MinTotCause + 1 
        # Generate the block cause matrix
        BlockCause[lower.tri(BlockCause)] = rnorm(n = nentries, mean = MeanCause, sd = SigmaCause);
        BlockCause[abs(BlockCause)>1] = 1;	# Don't use weights that are too high. 
      }
      Cause[AnchorInd[block]:(AnchorInd[block]+NModsInBlock[block]-1),
            AnchorInd[block]:(AnchorInd[block]+NModsInBlock[block]-1)] = BlockCause;
    }
    
    # Sprinkle a few more cross-block dependencies in
    nentries = sum(lower.tri(Cause))
    #NewEntries = runif(n = nentries, min = -MeanCause/2, max =  MeanCause/2);
    NewEntries = rnorm(n = nentries, sd =  MeanCause/2);
    # only keep about 1/2 of the new entries...
    #NewEntries[sample(nentries, size = as.integer(nentries/2))] = 0;
    Cause[lower.tri(Cause)] = Cause[lower.tri(Cause)] + NewEntries;
    # make sure anchors are independent
    Cause[AnchorInd, ] = 0;
    
    diag(Cause) = 0;
    
    NoiseVec = rep(Noise, times = No.Nodes); 
    
    # Generate the seed module eigengenes
    
    print.flush("Generating module eigengene seeds...");
    Seeds = CreateSeedVectors(Cause = Cause, AnchorInd = AnchorInd, AnchorVecs = AnchorVectors, 
                              Noise = NoiseVec);
    
    # Add a "grey seed vector" for compatibility with standard module finding stuff that assumes a grey
    # module is present
    
    ExtSeedVectors = cbind(Seeds$Vectors, rnorm(No.Samples));
    
    dev.set(3);
    image.plot(cor(Seeds$Vectors), col = GreenWhiteRed(50), zlim=c(-1,1))
    dev.set(2);
    
    # Generate the proportion of genes in each module
    
    Leftover = 1 - No.Nodes*MinProportion - GreyProportion;
    
    if (Leftover - GreyProportion<0.2) 
      stop("Not enough total genes for the number of modules and minimum simulated module size.");
    
    
    ModProps = rep(0, times = No.Nodes + 1);
    
    for (mod in 1:No.Nodes)
    {
      ModProps[mod] = MinProportion + 0.35*Leftover;
      Leftover = Leftover*0.55	# leave something for near-module genes
    }
    
    ModProps[No.Nodes + 1] = GreyProportion;
    
    # Mix up the module sizes so we don't get all the large modules in one block etc.
    
    ModProps[1:No.Nodes] = ModProps[sample(No.Nodes)];
    
    # Choose the modules to be left out in each set
    # LeaveOut contains entries for grey (always FALSE) for future use.
    
    LeaveOut = matrix(FALSE, nrow = No.Nodes, ncol = No.Sets);
    for (set in 1:No.Sets)
    {
      indices = sample(No.Nodes, size = as.integer(No.Nodes/5+1));
      LeaveOut[indices, set] = TRUE;
    }
    
    # Print some diagnostics
    
    print.flush(paste("SimulationLevel:", SimulationLevel, "simulation:", simulation));
    print.flush(paste("   Simulating", No.Nodes, " proper modules in", NModuleBlocks, 
                      "module blocks"));
    print.flush(paste("   Number of modules in each block: ", paste(NModsInBlock, collapse = ", ")));
    #print.flush(paste("   Number of genes in proper + grey modules: ", 
    #    paste(as.integer(ModProps * No.Genes), collapse = ", ")));
    print.flush(paste("   Adding", NOrderedSmallLayers[SimulationLevel],
                      " ordered and", NUnorderedSmallLayers[SimulationLevel],
                      " unordered layers, Random noise level:", BackgrNoise[SimulationLevel]));
    
    # Simulate the expression data.
    
    MinCorr = 0.6;
    NMods = No.Nodes;
    
    SimulatedData = SimulateSets(Seeds$Vectors, No.Genes, Subsets, ModProps, MinCorr, 
                                 BackgrNoise = BackgrNoise[SimulationLevel], 
                                 LeaveOut = LeaveOut,
                                 NOrderedSmallLayers = NOrderedSmallLayers[SimulationLevel] , 
                                 NUnorderedSmallLayers = NUnorderedSmallLayers[SimulationLevel] , 
                                 AverageNGenesInSmallModule = AverageNGenesInSmallModule[SimulationLevel] , 
                                 AverageExprInSmallModule = AverageExprInSmallModule[SimulationLevel] , 
                                 InvDensityOfSmallModules = 2, verbose = 2);
    
    ExprData = SimulatedData$ExprData;
    
    TrueColors = SimulatedData$TrueColors;
    TrueAllColors = SimulatedData$TrueAllColors;
    SimColorOrder = SimulatedData$ColorOrder;
    
    rm(SimulatedData);
    collect_garbage();
    
    #-------------------------------------------------------------------------------
    #
    # here comes the usual consensus analysis...
    #
    
    SoftPower = 6;
    
    KeepOverlapSign = FALSE;
    
    SelectedGenes = rep(TRUE, No.Genes);
    
    Connectivity = GetConnectivity(ExprData, Subsets, SoftPower, verbose = 1);
    collect_garbage();
    
    Modules = GetModules(ExprData, Subsets, SelectedGenes, SoftPower, 
                         BranchHeightCutoff = 0.95, ModuleMinSize = 30, 
                         KeepOverlapSign = KeepOverlapSign, verbose = 2);
    collect_garbage();
    
    Network = list(Subsets = Subsets, SoftPower = SoftPower, #Connectivity = Connectivity, 
                   SelectedGenes = SelectedGenes, Dissimilarity = Modules$Dissimilarity,
                   DissimilarityLevel = Modules$DissimilarityLevel, 
                   ClusterTree = Modules$ClusterTree, Colors = Modules$Colors,
                   BranchHeightCutoff = Modules$BranchHeightCutoff, 
                   ModuleMinSize = Modules$ModuleMinSize); 
    
    rm(Modules); collect_garbage();
    
    
    #-----------------------------------------------------------------------------------------------
    # Get intersection consensus modules
    
    ConsModMinSize = as.integer(MinSimModSize/6);
    
    Consensus = IntersectModules(Network = Network, 
                                 ConsBranchHeightCut = 0.97, ConsModMinSize = ConsModMinSize,
                                 verbose = 4)
    
    collect_garbage();
    DetectedColor = array(dim = c(No.Genes, No.DetBranchCuts, No.Sets+1));
    ColorOrder = order(SimColorOrder[, 1])
    ExtLeaveOut = rbind(LeaveOut, rep(FALSE, times = No.Sets));
    
    PlotCut = 3;
    cut = PlotCut;
    for (cut in 1:No.DetBranchCuts)
    {
      DetectedColor[ ,cut, 1] = as.character(cutreeDynamic(hierclust = Consensus$ClustTree, 
                                                           deepSplit = FALSE,
                                                           maxTreeHeight = DetBranchCuts[cut], minModuleSize=ConsModMinSize));
      for (set in 1:No.Sets)
      {
        DetectedColor[, cut, set+1] = as.character(cutreeDynamic(hierclust = Network$ClusterTree[[set]]$data, 
                                                                 deepSplit = FALSE,
                                                                 maxTreeHeight = DetBranchCuts[cut],
                                                                 minModuleSize = MinSimModSize/6));
        Network$Colors[, set] = DetectedColor[, cut, set+1];
      }
      # Attempt to assign simulated and consensus modules to one another.
      
      FTrueColors = as.factor(TrueAllColors[,1]); 
      No.TrueMods = nlevels(FTrueColors);
      for (set in 1:(No.Sets+1))
      {
        DetModColors = as.factor(DetectedColor[,cut,set]);
        No.DetMods = nlevels(DetModColors);
        pTable = matrix(0, nrow = No.TrueMods, ncol = No.DetMods);
        CountTbl = matrix(0, nrow = No.TrueMods, ncol = No.DetMods);
        for (smod in 1:No.TrueMods)
        {
          TrueMembers = (FTrueColors == levels(FTrueColors)[smod]);
          for (cmod in 1:No.DetMods)
          {
            DetMembers = (DetModColors == levels(DetModColors)[cmod]);
            if (No.DetMods>1)
            {
              pTable[smod, cmod] = fisher.test(TrueMembers, DetMembers, alternative = "greater")$p.value;
            } else {
              pTable[smod, cmod] = 1;
              if ((set==1) & (cut==PlotCut) & (PlotCut<No.DetBranchCuts)) PlotCut = PlotCut + 1;  
            }
            CountTbl[smod, cmod] = sum(FTrueColors == levels(FTrueColors)[smod] & 
                                         DetModColors == levels(DetModColors)[cmod])
          }
        }
        
        # Calculate module eigengenes for the found modules
        
        if (set==1)
        {
          PCs = NetworkModulePCs(ExprData, Network, UniversalModuleColors = DetectedColor[, cut, set],
                                 verbose=1)
        } else if (set==2) {
          PCs = NetworkModulePCs(ExprData, Network, verbose=1)
        }
        
        # Get the correlation of the seed vectors with the PCs
        
        if (set==1)
        {  
          SeedDetPCCor = cor(ExtSeedVectors[Network$Subsets==1, ColorOrder], PCs[[1]]$data);
          if (No.Sets>1) for (set2 in 2:No.Sets)
          {
            SeedDetPCCor = SeedDetPCCor + 
              cor(ExtSeedVectors[Network$Subsets==set2, ColorOrder], PCs[[set2]]$data);
          }
          SeedDetPCCor = SeedDetPCCor/No.Sets;
        } else {
          SeedDetPCCor = cor(ExtSeedVectors[Network$Subsets==set-1, ColorOrder], PCs[[set-1]]$data);
        }
        
        # Update the test diagnostics and statistics
        
        for (smod in 1:No.TrueMods)
        {
          ClosestCMod = order(pTable[smod, ])[1];
          NextClosest = order(pTable[smod, ])[2];
          count = sum(CountTbl[smod, ]);
          if (set==1)
          {
            SimStats$SimModIsConsensus[smod, set, simulation, SimulationLevel] = 
              (sum(ExtLeaveOut[ColorOrder[smod], ])==0) & (levels(FTrueColors)[smod]!="grey");
            if (SimStats$SimModIsConsensus[smod, set, simulation, SimulationLevel])
            {
              for (set2 in 1:No.Sets)
              {
                SModConn = Connectivity[TrueAllColors[,1]==levels(FTrueColors)[smod], set2];
                SModAssdConn = Connectivity[FTrueColors==levels(FTrueColors)[smod] &
                                              DetectedColor[,cut,1]!="grey", set2]; 
                SModUnassdConn = Connectivity[FTrueColors==levels(FTrueColors)[smod] &
                                                DetectedColor[,cut,1]=="grey", set2]; 
                SimStats$ConnMax[smod, set2, cut, simulation, SimulationLevel] = max(SModConn);
                SimStats$ConnAssdAve[smod, set2, cut, simulation, SimulationLevel] =
                  mean(SModAssdConn);
                SimStats$ConnUnassdAve[smod, set2, cut, simulation, SimulationLevel] =
                  mean(SModUnassdConn);
              }
            }
          } else {
            SimStats$SimModIsConsensus[smod, set, simulation, SimulationLevel] = 
              (ExtLeaveOut[ColorOrder[smod], set-1]==0) & (levels(FTrueColors)[smod]!="grey");
          }
          SimStats$SimModPValues[smod, set, cut, simulation, SimulationLevel] = pTable[smod, ClosestCMod]; 
          SimStats$SimModNextPValues[smod, set, cut, simulation, SimulationLevel] = pTable[smod, NextClosest];
          SimStats$HighestPValIsGrey[smod, set, cut, simulation, SimulationLevel] =
            (levels(DetModColors)[ClosestCMod]=="grey");
          SimStats$SimModPropInClosest[smod, set, cut, simulation, SimulationLevel] = 
            CountTbl[smod, ClosestCMod]/count;
          SimStats$SimModPropInGrey[smod, set, cut, simulation, SimulationLevel] = 
            CountTbl[smod, levels(DetModColors)=="grey"]/count;
          SimStats$SimModPropNotInGreyOrClosest[smod, set, cut, simulation, SimulationLevel] = 
            sum(CountTbl[smod, levels(DetModColors)!="grey" &
                           levels(DetModColors)!=levels(DetModColors)[ClosestCMod]])/count;
          SimStats$Seed_ClosestPC_Cor[smod, set, cut, simulation, SimulationLevel] = 
            SeedDetPCCor[smod, ClosestCMod];
        }
        if (cut==PlotCut & set==1) 
        {  
          PlotpTable = pTable; 
          PlotCountTable = CountTbl; 
          PlotSeedPCCor = SeedDetPCCor;
          PlotDetModColors = DetModColors;
        }
        if (set==1)
        {
          for (set2 in 1:No.Sets)
          {
            SimStats$PCCor[set2, cut, simulation, SimulationLevel] = list(data = cor(PCs[[set]]$data));
          }
        }
      }
    }
    
    collect_garbage();
    if (simulation <= PlotsPerLevel)
    {
      # Redo cutting of the set clustering trees by dynamic cutting with fewer branch cuts.
      
      ModuleMinSize = as.integer(MinSimModSize/5*3);
      SetColors = array("grey", dim = c(No.Genes, No.SetBranchCuts*No.Sets));
      SetPlotLabels = NULL;
      for (set2 in 1:No.Sets)
      {
        for (setcut in 1:No.SetBranchCuts)
        {
          SetColors[, (set2-1)*No.SetBranchCuts + setcut] = 
            as.character(cutreeDynamic(hierclust = Network$ClusterTree[[set2]]$data, 
                                       deepSplit = FALSE,
                                       maxTreeHeight = SetBranchCuts[setcut], minModuleSize=ModuleMinSize));
          SetPlotLabels = c(SetPlotLabels, paste("Set", set2, "@", SetBranchCuts[setcut]));
        }
      }
      
      collect_garbage();
      
      # Plot the consensus dendrogram and the corresponding colors
      
      StandardCex = 1;
      FileNumberID = paste("Lev-", SimulationLevel, "-No-", simulation, "-", sep="");
      SizeWindow(12,10);
      par(mfrow = c(2,1));
      par(mar=c(3,5,2,1.2)+0.2);
      
      plot(Consensus$ClustTree, labels = FALSE, main = "Consensus dendrogram");
      abline(DetBranchCuts[PlotCut], 0, col="red");
      hclustplotn(Consensus$ClustTree, cbind(TrueAllColors[,1], SetColors, DetectedColor[, , 1]), 
                  RowLabels = c("Sim\'d", SetPlotLabels, paste("Cons", DetBranchCuts)),
                  main=paste("Module colors"));
      
      dev.copy2eps(file = paste(PlotDir, PlotFileBase, FileNumberID, "DetDendr.eps", sep = ""));
      
      # Plot the new module membership vs. set dendrograms
      
      SizeWindow(12,10);
      par(mar=c(5,5,2,2)+0.1);
      par(mfcol=c(2,2));
      for (i in (1:No.Sets))
      {
        plot(Network$ClusterTree[[i]]$data,labels=F,xlab="",main=paste("Dendrogram", Set.Labels[i]),
             ylim=c(0,1), sub = "")
        SCColor = cbind(TrueAllColors[, 1], SetColors, DetectedColor); 
        RowLabels = c("Sim\'d", SetPlotLabels, paste("Cons", DetBranchCuts));
        hclustplotn(Network$ClusterTree[[i]]$data, SCColor, RowLabels = RowLabels, cex.RowLabels = 1,
                    main="Module colors")
      }
      
      dev.copy2eps(file = paste(PlotDir, PlotFileBase, FileNumberID, "SetDendr.eps", sep = ""));
      
      # Plot assignment p-values and correlation of PCs with seeds
      
      SizeWindow(12,6);
      par(mfrow=c(1,2));
      par(cex = StandardCex);
      par(mar=c(3,5,3,3)+0.3);
      
      LogpTable = -log(PlotpTable);
      LogpTable[is.na(LogpTable)] = 50;
      LogpTable[LogpTable>50] = 50;
      SimSymbol = ifelse(SimStats$SimModIsConsensus[, 1, simulation, SimulationLevel], "Con", "XXX");
      HeatmapWithTextLabels(Matrix = LogpTable,
                            xLabels = paste("XX", levels(PlotDetModColors), sep=""),
                            yLabels = paste("XX", levels(FTrueColors), sep=""),
                            ySymbols = SimSymbol,
                            NumMatrix = PlotCountTable, ColorLabels = TRUE,
                            InvertColors = TRUE, SetMargins = FALSE,
                            main = "-log(p) (by color) and counts");
      
      HeatmapWithTextLabels(Matrix = PlotSeedPCCor,
                            colors = GreenWhiteRed(50), zlim = c(-1,1),
                            xLabels = paste("XX", levels(PlotDetModColors), sep=""),
                            yLabels = paste("XX", levels(FTrueColors), sep=""),
                            ySymbols = SimSymbol,
                            NumMatrix = signif(PlotSeedPCCor*100,2), ColorLabels = TRUE,
                            InvertColors = FALSE, SetMargins = FALSE,
                            main = "Correlation (%)");
      
      dev.copy2eps(file = paste(PlotDir, PlotFileBase, FileNumberID, "SimVsCalcMods.eps", sep = ""));
    }
    collect_garbage();
  }
}

save(SimStats, file = paste(PlotFileBase, "SimStats.RData", sep=""));

# Stop here

stop("Simulation is now done.");

# Process simulation data

#load(file = paste(PlotFileBase, "SimStats.RData", sep=""));

NoMods = array(0, dim = c(2, 2, No.Sets+1, No.DetBranchCuts, No.SimulationLevels));

for (level in 1:No.SimulationLevels)
{
  for (cut in 1:No.DetBranchCuts)
  {
    for (set in 1:(No.Sets+1))
    {
      for (SimMod in c(0,1)) for (DetMod in c(0,1))
      {
        NoMods[SimMod+1, DetMod+1, set, cut, level] = 
          sum(SimStats$SimModIsConsensus[ , set, , level]==SimMod & 
                SimStats$HighestPValIsGrey[ , set, cut, , level]==(1-DetMod), na.rm = TRUE);
      } 
      print(paste("Level: ", level, ", cut: ", DetBranchCuts[cut], ", set: ", set-1, 
                  ", S&D: ", NoMods[2, 2, set, cut, level], ", S&ND: ", NoMods[2,1,set,cut,level],
                  ", NS&D: ", NoMods[1, 2, set, cut, level], ", NS&ND: ", NoMods[1,1,set,cut,level], 
                  sep = ""));
    }
  }
}

# Write the same information in a latex table?

LatexFile = paste(PlotFileBase, "SensitivityTable.tex", sep="");

# These are the table entries
NoiseLevel = NULL; Cut = NULL; ConsSens = NULL; ConsSpec = NULL; SetSens = NULL; SetSpec = NULL;
ConnAssd = NULL; ConnUnassd = NULL;

for (level in 1:No.SimulationLevels)
{
  for (cut in 1:No.DetBranchCuts)
  {
    NoiseLevel = c(NoiseLevel, level);
    Cut = c(Cut, DetBranchCuts[cut]);
    ConsSens = c(ConsSens, signif(NoMods[2, 2, 1, cut, level]/(NoMods[2, 2, 1, cut, level] + 
                                                                 NoMods[2, 1, 1, cut, level]),3));
    ConsSpec = c(ConsSpec, signif(NoMods[2, 2, 1, cut, level]/(NoMods[2, 2, 1, cut, level] + 
                                                                 NoMods[1, 2, 1, cut, level]),3));
    SetSD = sum(NoMods[2, 2, c(2:(No.Sets+1)), cut, level]);
    SetNSD = sum(NoMods[1, 2, c(2:(No.Sets+1)), cut, level]);
    SetSND = sum(NoMods[2, 1, c(2:(No.Sets+1)), cut, level]);
    SetSens = c(SetSens, signif(SetSD/(SetSD + SetSND),3));
    SetSpec = c(SetSpec, signif(SetSD/(SetSD + SetNSD),3));
  }
}

SimModIsConsensus_Ext = SimStats$HighestPValIsGrey; # To get the dimensions right
for (cut in 1:No.DetBranchCuts)
  SimModIsConsensus_Ext[, , cut, , ] = SimStats$SimModIsConsensus[,,,];

Quantile = 0.95;
#ConsRatioAboveQuantile = matrix(0, nrow= No.DetBranchCuts, ncol = No.SimulationLevels);
#SetRatioAboveQuantile = matrix(0, nrow= No.DetBranchCuts, ncol = No.SimulationLevels);
ConsRatioAboveQuantile = NULL;
SetRatioAboveQuantile = NULL;
for (level in 1:No.SimulationLevels)
{
  for (cut in 1:No.DetBranchCuts)
  {
    Include = SimModIsConsensus_Ext[,1,cut,,level] & (SimStats$HighestPValIsGrey[,1,cut,,level]==FALSE);
    Seed_PCCor = SimStats$Seed_ClosestPC_Cor[,1,cut,,level];
    ConsRatioAboveQuantile = c(ConsRatioAboveQuantile, 
                               sum(abs(Seed_PCCor[Include])>Quantile, na.rm = TRUE)/sum(Include==TRUE, na.rm = TRUE));
    SumAbove = 0; SumTotal = 0;
    for (set in 2:(No.Sets+1))
    {
      Include = SimModIsConsensus_Ext[,set,cut,,level] & (SimStats$HighestPValIsGrey[,set,cut,,level]==FALSE);
      Seed_PCCor = SimStats$Seed_ClosestPC_Cor[,set,cut,,level];
      SumAbove = SumAbove + sum(abs(Seed_PCCor[Include])>Quantile, na.rm = TRUE)
      SumTotal = SumTotal + sum(Include==TRUE, na.rm = TRUE);
    }
    SetRatioAboveQuantile = c(SetRatioAboveQuantile, SumAbove/SumTotal);
  }
}

LatexTable = data.frame(NoiseLevel, Cut, ConsSens, ConsSpec, signif(ConsRatioAboveQuantile,3),
                        SetSens, SetSpec, signif(SetRatioAboveQuantile,3)); 
names(LatexTable) = c("NL", "Cut", "Sensitivity", "Specificity", "$P_{0.95}$", 
                      "Sensitivity", "Specificity",  "$P_{0.95}$");
write.table(LatexTable, file=LatexFile, row.names=FALSE, col.names=TRUE, quote = FALSE, sep = " & ", 
            eol = "\\\\ \n ");

ConnAssd = SimStats$ConnAssdAve;
ConnAssd[is.na(ConnAssd)] = 0;
ConnUnassd = SimStats$ConnUnassdAve;
ConnUnassd[is.na(ConnUnassd)] = 0;
mx = max(ConnAssd);
SizeWindow(12,9)
par(mfrow=c(2,3))
par(cex = 1.0);
for (level in 1:No.SimulationLevels)
{
  plot(as.vector(SimStats$ConnAssdAve[,,,,level]), as.vector(SimStats$ConnUnassdAve[,,,,level]), 
       # xlim = c(0,mx), ylim = c(0,mx), 
       cex = 0.6, cex.main = 1.2, cex.axis = 1.0, cex.lab = 1.0,
       xlab = "K.in", 
       ylab = "K.out",
       main = paste("Noise level:", level));
  abline(0,1,col = "blue");
}

dev.copy2eps(file = paste(PlotFileBase, "AssignedVsUnassignedConn.eps", sep = ""));


SizeWindow(12,9)
par(mfrow=c(2,3))
par(cex = 1.0);
for (level in 1:No.SimulationLevels)
{
  print(level)
  IncludeInHistogram = SimModIsConsensus_Ext[,1,,,level] & (SimStats$HighestPValIsGrey[,1,,,level]==FALSE);
  Seed_PCCor = SimStats$Seed_ClosestPC_Cor[,1,,,level];
  # Seed_PCCor = as.vector(Seed_PCCor[as.vector(!is.na(Seed_PCCor))]);
  h = hist(abs(Seed_PCCor[IncludeInHistogram]), main = paste("Noise level:", level), 
           xlab = "Correlation", breaks = 400 );
}
dev.copy2eps(file = paste(PlotFileBase, "Seed-ClosestPC-Cor.eps", sep = ""));

