loadModule("som_module", TRUE)
RcppParallel::setThreadOptions(numThreads = RcppParallel::defaultNumThreads() - 1)
#loadModule("trg_module", TRUE)