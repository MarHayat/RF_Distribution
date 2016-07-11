require(ape)
require(igraph)

MaxLOG = 700

driverBeast = function(Collection, complexTable, baseFilename, initTree = NULL, analysis = FALSE, repeated = FALSE) {
  treeFilename = paste(baseFilename, ".trees", sep = "")
  if (analysis) {
    output = processBeastOutput(treeFilename, Collection, complexTable, analysis = TRUE)
  }
  else if (repeated) {
    output = processBeastOutput(treeFilename, Collection, complexTable, analysis = FALSE)
  }
  else {
    Input = createBeastInput(Collection)
    writeBeastXML(Input, baseFilename, initTree = initTree)
    fullFilename = paste(baseFilename, ".xml", sep = "")
    initDir = getwd()
    fullFilename = paste(initDir, "/", fullFilename, sep = "")
    setwd(beastDir)
    system(paste(c("./beast", "-threads", "0", fullFilename), collapse = " "))
    ### check to make sure the file is written before copying it!
    logFilename  = paste(baseFilename, ".log",   sep = "")
    while (!all(file.exists(c(treeFilename, logFilename)))) {
      print("...")
      U = runif(1e6)
    }
    file.copy(from = c(treeFilename, logFilename), to = initDir)
    file.remove(c(treeFilename, logFilename))
    setwd(initDir)
    output = processBeastOutput(treeFilename, Collection, complexTable)
    outFilename = paste(baseFilename, ".results", sep = "")
    createDescription(output[[2]], complexTable, outFilename)
    save(output, file = paste(baseFilename, ".RData", sep = ""))
  }
  output
}

processDataFiles = function(optionsList, baseline = "BeastMIRU_Latest", extension = ".RData", special = 3) {
  FUN = function(x, y) {outer(x, y, paste0)}
  L = length(optionsList)
  fullOptions = optionsList[[L]]
  for (ind in (L-1):1) {
    fullOptions = as.vector(do.call(FUN, list(optionsList[[ind]], fullOptions)))
  }
  Q = length(fullOptions)
  Env = new.env()
  for (opt in 1:Q) {
    curOpt = fullOptions[opt]
    print(curOpt)
    curFile = paste0(baseline, curOpt, extension)
    load(curFile, Env)
    curPreds = get("output", Env)[[2]]
    if (opt == 1) {
      Results = matrix(NA, length(curPreds), Q)
      rownames(Results) = names(curPreds)
      colnames(Results) = fullOptions
      Results[,1] = curPreds
    }
    else {
      Results[names(curPreds), opt] = curPreds	
    }
  }
  outFilename = paste0(baseline, optionsList[[special]][1], ".data")
  write.table(Results, file = outFilename, quote = FALSE, sep = "\t")
  Results
}

createDescription = function(predictions, complexTable, outFile) {
  L = length(complexTable)
  Names = names(predictions)
  Lines = rep(NA, L)
  Lengths = sapply(complexTable, length)
  for (ind in 1:L) {
    Dim = Lengths[ind]
    Prob = 1 - predictions[ind]
    Name = Names[ind]
    Lines[ind] = paste("Patient", Name, "has", Dim, "strains and a mixed infection with probability", Prob)
  }
  f = file(outFile, "w")
  writeLines(Lines, f)
  close(f)
}

createBeastInput = function(Collection, base = "C", other = "T") {
  m = nrow(Collection)
  n = ncol(Collection)
  N = n * (maxCNV + 1)
  Input = matrix(base, m, N)
  for (ind1 in 1:m) {
    for (ind2 in 1:n) {
      curEntry = Collection[ind1, ind2]
      if (curEntry > 0) {
        startInd  = (ind2 - 1) * (maxCNV + 1) + 1
        finishInd = startInd - 1 + curEntry
        curRange = startInd:finishInd
        Input[ind1, curRange] = other
      }
    }
  }
  rownames(Input) = rownames(Collection)
  Input
}

writeBeastFile = function(Matrix, outFile) {
  m = nrow(Matrix)
  n = ncol(Matrix)
  Lines = rep('', m + 7)
  Lines[1] = "#NEXUS"
  Lines[2] = "Begin data;"
  Lines[3] = paste("Dimensions ntax=", m," nchar=", n, ";", sep = "")
  Lines[4] = "Format datatype=dna missing=? gap=-;"
  Lines[5] = "Matrix"
  for (ind in 1:m) {
    Lines[5 + ind] = paste(rownames(Matrix)[ind], paste(Matrix[ind,], collapse = ""), sep = "\t")
  }
  Lines[m + 6] = ";"
  Lines[m + 7] = "End"
  f = file(outFile, "w")
  writeLines(Lines, f)
  close(f)
}

writeBeastXML = function(Matrix, baseFilename, initTree, templateFilename = "BeastProcessing") {
  m = nrow(Matrix)
  n = ncol(Matrix)
  Names = rownames(Matrix)
  Lines = rep('', 5 * m + 10)
  Lines[1] = '<?xml version="1.0" standalone="yes"?>'
  Lines[3] = '<beast>'
  Lines[5] = '\t<taxa id="taxa">'
  for (ind in 1:m) {
    Lines[5 + ind] = paste('\t\t<taxon id="', Names[ind],'"/>', sep = "")
  }
  Lines[m + 6] = '\t</taxa>'
  Lines[m + 8] = '\t<alignment id="alignment" dataType="nucleotide">'
  for (ind in 1:m) {
    Lines[m + 8 + 4 * (ind - 1) + 1] = '\t\t<sequence>'
    Lines[m + 8 + 4 * (ind - 1) + 2] = paste('\t\t\t<taxon idref="', Names[ind], '"/>', sep = "")
    Lines[m + 8 + 4 * (ind - 1) + 3] = paste('\t\t\t\t', paste(Matrix[ind,], collapse = ""), sep = "")
    Lines[m + 8 + 4 * (ind - 1) + 4] = '\t\t</sequence>'
  }
  Lines[5 * m + 9] = '\t</alignment>'
  auxFile = paste(templateFilename, ".xml", sep = "")
  f = file(auxFile, "r")
  otherLines = readLines(f)
  close(f)
  otherLines = gsub(templateFilename, baseFilename, otherLines, fixed = TRUE) # replace filenames as needed!
  if (!is.null(initTree)) {
    phyloTree = convertPhylo(initTree)
    treeFile = paste(baseFilename, "InitTree", ".nex", sep = "")
    write.nexus(phyloTree, file = treeFile)
    fr = file(treeFile, "r")
    lines = readLines(fr)
    close(fr)
    goodLine = which(nchar(lines) > nchar(paste(Names, collapse = "")))
    treeLine = lines[goodLine]
    splitLine = unlist(strsplit(treeLine, " "))
    pureTree = paste('\t\t', splitLine[length(splitLine)], sep = "")
    substLine = grep('orangutan', otherLines)
    if (length(substLine) != 1) {
      stop("Error: can't find a line to replace with the user-defined tree")
    }
    otherLines[substLine] = pureTree
  }
  else { # eventually make provisions for using UPGMA tree
    print("Warning: no starting tree has been defined!")
  }
  outFile = paste(baseFilename, ".xml", sep = "")
  allLines = c(Lines, otherLines)
  fw = file(outFile, "w")
  writeLines(allLines, fw)
  close(fw)
}

processBeastOutput = function(filename, collection, complexTable, burnIn = 0.1, analysis = FALSE) {
  f = file(filename, "r")
  allLines = readLines(f)
  close(f)
  if (length(allLines) == 0) {
    print("Warning: the file is empty!")
    return(list(NULL, NULL))
  }
  Probs = regexpr("posterior=-[0-9.]*", allLines)
  Matches = substr(allLines, as.integer(Probs), as.integer(Probs) + attr(Probs, "match.length") - 1)
  Matches = Matches[nchar(Matches) > 0]
  Probs = substr(Matches, nchar("posterior=-"), nchar(Matches))
  Probs = as.numeric(Probs)
  L0 = length(Probs)
  allTrees = read.nexus(filename)
  if (length(allTrees) != L0) {
    stop("Error: size mismatch between trees and probabilities!")
  }
  if (burnIn > 0) {
    skipTrees = round(burnIn * L0)
    Probs = Probs[-(1:skipTrees)]
    allTrees = allTrees[-(1:skipTrees)]
    print(paste("Removing", skipTrees, "out of", L0, "to account for burn-in"))
    L0 = L0 - skipTrees 
  }
  Probs = Probs - max(Probs) # shift probabilities to have max 0
  goodInds = (Probs > - MaxLOG)
  redProbs = Probs[goodInds]
  redTrees = allTrees[goodInds]
  L = length(redProbs)
  print(paste("Only", L, "out of", L0, "trees have non-negligible probability"))
  numComplex = length(complexTable)
  if (analysis) {
    output = vector("list", L)
    names(output) = redProbs
  }
  else {
    ProbMatrix = matrix(NA, L, numComplex)
    rownames(ProbMatrix) = redProbs
    colnames(ProbMatrix) = names(complexTable)
  }
  for (index in 1:L) {
    if (index %% 100 == 0) {
      print(paste("Processed", index, "trees so far"))
    }
    curTree = convertTree(redTrees[[index]])
    curLeaves = (length(curTree) + 1)/2
    curLeafNames = as.numeric(names(curTree[1:curLeaves]))
    curCollection = collection[curLeafNames, ]
    reconstruction = intervalReconstruct(curCollection, curTree)
    if (analysis) {
      output[[index]] = vector("list", numComplex)
      names(output[[index]]) = names(complexTable)
    }
    for (ind in 1:numComplex) {
      indices = complexTable[[ind]]
      trueIndices = which(curLeafNames %in% indices)
      stopifnot(length(indices) == length(trueIndices))
      if (analysis) {
        output[[index]][[ind]] = analyzeTree(reconstruction, curTree, trueIndices)
      }
      else {
        ProbMatrix[index, ind] = sum(getMultiProbabilityEvolved(reconstruction, curTree, trueIndices))
      }
    }
  }
  if (!analysis) {
    print("These are the distinct values in the probability matrix")
    print(unique(as.vector(ProbMatrix)))
    if (any(ProbMatrix > 1)) {
      print("Warning: some probabilities exceed 1; this should never happen!")
    }
    predictions = colSums(exp(redProbs) * ProbMatrix) / sum(exp(redProbs))
    output = list(ProbMatrix, predictions)
  }
  output
}

checkFile = function(filename, burnIn = 0.1) {
  f = file(filename, "r")
  allLines = readLines(f)
  close(f)
  Probs = regexpr("posterior=-[0-9.]*", allLines)
  ProbsQ = regexpr("posterior=-[0-9.]*E[0-9.]", allLines)
  baseMatches = which(Probs >= 0)
  QMatches = which(ProbsQ >= 0)
  L = length(baseMatches)
  cutoff = baseMatches[as.integer(L * burnIn)]
  if (length(QMatches) > 0 && max(QMatches) > cutoff) {
    return(FALSE)
  }
  else {
    return(TRUE)
  }
}

checkFiles = function(optionsList, baseline = "BeastMIRU_Latest", extension = ".trees") {
  FUN = function(x, y) {outer(x, y, paste0)}
  L = length(optionsList)
  fullOptions = optionsList[[L]]
  for (ind in (L-1):1) {
    fullOptions = as.vector(do.call(FUN, list(optionsList[[ind]], fullOptions)))
  }
  Q = length(fullOptions)
  checks = rep(NA, Q)
  names(checks) = fullOptions
  for (opt in 1:Q) {
    curOpt = fullOptions[opt]
    print(curOpt)
    curFile = paste0(baseline, curOpt, extension)
    checks[opt] = checkFile(curFile)
  }
  if (!all(checks)) {
    print(checks)
    stop("Error: not all checks have been passed!")
  }
}