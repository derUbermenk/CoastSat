import sys

module = type(sys)("coastsat")
sys.modules["coastsat"] = module

module.submodule = type(sys)("SDS_download")
sys.modules["coastsat.SDS_download"] = module.submodule

module.submodule = type(sys)("SDS_shoreline")
sys.modules["coastsat.SDS_shoreline"] = module.submodule

module.submodule = type(sys)("SDS_tools")
sys.modules["coastsat.SDS_tools"] = module.submodule

module.submodule = type(sys)("SDS_transects")
sys.modules["coastsat.SDS_transects"] = module.submodule