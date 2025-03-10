#include "DataFormats/Provenance/interface/ProductDescription.h"

#include "FWCore/Utilities/interface/EDMException.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "FWCore/Utilities/interface/FriendlyName.h"
#include "FWCore/Reflection/interface/FunctionWithDict.h"
#include "FWCore/Reflection/interface/TypeWithDict.h"
#include "FWCore/Utilities/interface/WrappedClassName.h"

#include "TDictAttributeMap.h"

#include <cassert>
#include <ostream>
#include <sstream>

class TClass;

namespace edm {
  BranchDescription::Transients::Transients()
      : branchName_(),
        wrappedName_(),
        wrappedType_(),
        unwrappedType_(),
        produced_(false),
        onDemand_(false),
        isTransform_(false),
        dropped_(false),
        transient_(false),
        availableOnlyAtEndTransition_(false),
        isMergeable_(false) {}

  void BranchDescription::Transients::reset() { *this = BranchDescription::Transients(); }

  BranchDescription::BranchDescription()
      : branchType_(InEvent),
        moduleLabel_(),
        processName_(),
        branchID_(),
        fullClassName_(),
        friendlyClassName_(),
        productInstanceName_(),
        branchAliases_(),
        aliasForBranchID_(),
        transient_() {
    // do not call init here! It will result in an exception throw.
  }

  BranchDescription::BranchDescription(BranchType const& branchType,
                                       std::string const& moduleLabel,
                                       std::string const& processName,
                                       std::string const& className,
                                       std::string const& friendlyClassName,
                                       std::string const& productInstanceName,
                                       TypeWithDict const& theTypeWithDict,
                                       bool produced,
                                       bool availableOnlyAtEndTransition,
                                       std::set<std::string> const& aliases)
      : branchType_(branchType),
        moduleLabel_(moduleLabel),
        processName_(processName),
        branchID_(),
        fullClassName_(className),
        friendlyClassName_(friendlyClassName),
        productInstanceName_(productInstanceName),
        branchAliases_(aliases),
        transient_() {
    setDropped(false);
    setProduced(produced);
    setOnDemand(false);
    transient_.availableOnlyAtEndTransition_ = availableOnlyAtEndTransition;
    setUnwrappedType(theTypeWithDict);
    init();
  }

  BranchDescription::BranchDescription(BranchDescription const& aliasForBranch,
                                       std::string const& moduleLabelAlias,
                                       std::string const& productInstanceAlias)
      : branchType_(aliasForBranch.branchType()),
        moduleLabel_(moduleLabelAlias),
        processName_(aliasForBranch.processName()),
        branchID_(),
        fullClassName_(aliasForBranch.className()),
        friendlyClassName_(aliasForBranch.friendlyClassName()),
        productInstanceName_(productInstanceAlias),
        branchAliases_(aliasForBranch.branchAliases()),
        aliasForBranchID_(aliasForBranch.branchID()),
        transient_() {
    setDropped(false);
    setProduced(aliasForBranch.produced());
    setOnDemand(false);  // will be re-set externally to the aliasForBranch.onDemand() after that one has been set
    transient_.availableOnlyAtEndTransition_ = aliasForBranch.availableOnlyAtEndTransition();
    setUnwrappedType(aliasForBranch.unwrappedType());
    init();
  }

  void BranchDescription::initBranchName() {
    if (!branchName().empty()) {
      return;  // already called
    }
    throwIfInvalid_();

    char const underscore('_');
    char const period('.');

    if (friendlyClassName_.find(underscore) != std::string::npos) {
      throw cms::Exception("IllegalCharacter")
          << "Class name '" << friendlyClassName()
          << "' contains an underscore ('_'), which is illegal in the name of a product.\n";
    }

    // Module labels of non-persistent products are allowed to contain
    // underscores. For module labels of persistent products, the module
    // label is checked for underscores in the function initFromDictionary
    // after we determine whether the product is persistent or not.

    if (productInstanceName_.find(underscore) != std::string::npos) {
      throw cms::Exception("IllegalCharacter")
          << "Product instance name '" << productInstanceName()
          << "' contains an underscore ('_'), which is illegal in a product instance name.\n";
    }

    if (processName_.find(underscore) != std::string::npos) {
      throw cms::Exception("IllegalCharacter")
          << "Process name '" << processName()
          << "' contains an underscore ('_'), which is illegal in a process name.\n";
    }

    std::string& brName = transient_.branchName_;
    brName.reserve(friendlyClassName().size() + moduleLabel().size() + productInstanceName().size() +
                   processName().size() + 4);
    brName += friendlyClassName();
    brName += underscore;
    brName += moduleLabel();
    brName += underscore;
    brName += productInstanceName();
    brName += underscore;
    brName += processName();
    brName += period;

    if (!branchID_.isValid()) {
      branchID_.setID(brName);
    }
  }

  void BranchDescription::initFromDictionary() {
    if (bool(wrappedType())) {
      return;  // already initialized;
    }

    throwIfInvalid_();

    try {
      setWrappedName(wrappedClassName(fullClassName()));
      // unwrapped type.
      setUnwrappedType(TypeWithDict::byName(fullClassName()));
      if (!bool(unwrappedType())) {
        setTransient(false);
        return;
      }
    } catch (edm::Exception& caughtException) {
      caughtException.addContext(std::string{"While initializing meta data for branch: "} + branchName());
      throw;
    }

    edm::TypeWithDict wrType(TypeWithDict::byName(wrappedName()));
    try {
      setWrappedType(wrType);
      if (!bool(wrappedType())) {
        return;
      }
    } catch (edm::Exception& caughtException) {
      caughtException.addContext(std::string{"While initializing meta data for branch: "} + branchName());
      throw;
    }

    setTransient(false);
    TDictAttributeMap* wp = wrappedType().getClass()->GetAttributeMap();
    if (wp && wp->HasKey("persistent") && !strcmp(wp->GetPropertyAsString("persistent"), "false")) {
      // Set transient if persistent == "false".
      setTransient(true);
      return;
    } else {
      // Module labels of persistent products cannot contain underscores,
      // but for non-persistent products it is allowed because path names
      // are used as module labels for path status products and there
      // are many path names that include underscores.
      char const underscore('_');
      if (moduleLabel_.find(underscore) != std::string::npos) {
        throw cms::Exception("IllegalCharacter")
            << "Module label '" << moduleLabel()
            << "' contains an underscore ('_'), which is illegal in a module label.\n";
      }
    }
  }

  void BranchDescription::merge(BranchDescription const& other) {
    branchAliases_.insert(other.branchAliases().begin(), other.branchAliases().end());
  }

  void BranchDescription::setSwitchAliasForBranch(BranchDescription const& aliasForBranch) {
    if (branchType_ != aliasForBranch.branchType()) {
      throw Exception(errors::LogicError) << "BranchDescription::setSwitchAliasForBranch: branchType (" << branchType_
                                          << ") differs from aliasForBranch (" << aliasForBranch.branchType()
                                          << ").\nPlease report this error to the FWCore developers";
    }
    if (produced() != aliasForBranch.produced()) {
      throw Exception(errors::LogicError) << "BranchDescription::setSwitchAliasForBranch: produced differs from "
                                             "aliasForBranch.\nPlease report this error to the FWCore developers";
    }
    if (unwrappedTypeID().typeInfo() != aliasForBranch.unwrappedType().typeInfo()) {
      throw Exception(errors::LogicError)
          << "BranchDescription::setSwitchAliasForBranch: unwrapped type info (" << unwrappedTypeID().name()
          << ") differs from aliasForBranch (" << aliasForBranch.unwrappedType().typeInfo().name()
          << ").\nPlease report this error to the FWCore developers";
    }

    branchAliases_ = aliasForBranch.branchAliases();
    transient_.switchAliasForBranchID_ = aliasForBranch.originalBranchID();
    transient_.availableOnlyAtEndTransition_ = aliasForBranch.availableOnlyAtEndTransition();
  }

  void BranchDescription::write(std::ostream& os) const {
    os << "Branch Type = " << branchType() << std::endl;
    os << "Process Name = " << processName() << std::endl;
    os << "ModuleLabel = " << moduleLabel() << std::endl;
    os << "Branch ID = " << branchID() << '\n';
    os << "Class Name = " << fullClassName() << '\n';
    os << "Friendly Class Name = " << friendlyClassName() << '\n';
    os << "Product Instance Name = " << productInstanceName() << std::endl;
  }

  void throwExceptionWithText(char const* txt) {
    Exception e(errors::LogicError);
    e << "Problem using an incomplete BranchDescription\n"
      << txt << "\nPlease report this error to the FWCore developers";
    throw e;
  }

  void BranchDescription::throwIfInvalid_() const {
    if (branchType_ >= NumBranchTypes)
      throwExceptionWithText("Illegal BranchType detected");

    if (moduleLabel_.empty())
      throwExceptionWithText("Module label is not allowed to be empty");

    if (processName_.empty())
      throwExceptionWithText("Process name is not allowed to be empty");

    if (fullClassName_.empty())
      throwExceptionWithText("Full class name is not allowed to be empty");

    if (friendlyClassName_.empty())
      throwExceptionWithText("Friendly class name is not allowed to be empty");
  }

  void BranchDescription::updateFriendlyClassName() {
    friendlyClassName_ = friendlyname::friendlyName(fullClassName());
    clearBranchName();
    initBranchName();
  }

  bool operator<(BranchDescription const& a, BranchDescription const& b) {
    if (a.processName() < b.processName())
      return true;
    if (b.processName() < a.processName())
      return false;
    if (a.fullClassName() < b.fullClassName())
      return true;
    if (b.fullClassName() < a.fullClassName())
      return false;
    if (a.friendlyClassName() < b.friendlyClassName())
      return true;
    if (b.friendlyClassName() < a.friendlyClassName())
      return false;
    if (a.productInstanceName() < b.productInstanceName())
      return true;
    if (b.productInstanceName() < a.productInstanceName())
      return false;
    if (a.moduleLabel() < b.moduleLabel())
      return true;
    if (b.moduleLabel() < a.moduleLabel())
      return false;
    if (a.branchType() < b.branchType())
      return true;
    if (b.branchType() < a.branchType())
      return false;
    if (a.branchID() < b.branchID())
      return true;
    if (b.branchID() < a.branchID())
      return false;
    if (a.branchAliases() < b.branchAliases())
      return true;
    if (b.branchAliases() < a.branchAliases())
      return false;
    if (a.present() < b.present())
      return true;
    if (b.present() < a.present())
      return false;
    return false;
  }

  bool combinable(BranchDescription const& a, BranchDescription const& b) {
    return (a.branchType() == b.branchType()) && (a.processName() == b.processName()) &&
           (a.fullClassName() == b.fullClassName()) && (a.friendlyClassName() == b.friendlyClassName()) &&
           (a.productInstanceName() == b.productInstanceName()) && (a.moduleLabel() == b.moduleLabel()) &&
           (a.branchID() == b.branchID());
  }

  bool operator==(BranchDescription const& a, BranchDescription const& b) {
    return combinable(a, b) && (a.dropped() == b.dropped()) && (a.branchAliases() == b.branchAliases());
  }

  std::string match(BranchDescription const& a, BranchDescription const& b, std::string const& fileName) {
    std::ostringstream differences;
    if (a.branchName() != b.branchName()) {
      differences << "Branch name '" << b.branchName() << "' does not match '" << a.branchName() << "'.\n";
      // Need not compare components of branch name individually.
      // (a.friendlyClassName() != b.friendlyClassName())
      // (a.moduleLabel() != b.moduleLabel())
      // (a.productInstanceName() != b.productInstanceName())
      // (a.processName() != b.processName())
    }
    if (a.branchType() != b.branchType()) {
      differences << "Branch '" << b.branchName() << "' is a(n) '" << b.branchType() << "' branch\n";
      differences << "    in file '" << fileName << "', but a(n) '" << a.branchType()
                  << "' branch in previous files.\n";
    }
    if (a.branchID() != b.branchID()) {
      differences << "Branch '" << b.branchName() << "' has a branch ID of '" << b.branchID() << "'\n";
      differences << "    in file '" << fileName << "', but '" << a.branchID() << "' in previous files.\n";
    }
    if (a.fullClassName() != b.fullClassName()) {
      differences << "Products on branch '" << b.branchName() << "' have type '" << b.fullClassName() << "'\n";
      differences << "    in file '" << fileName << "', but '" << a.fullClassName() << "' in previous files.\n";
    }
    if (!b.dropped() && a.dropped()) {
      differences << "Branch '" << a.branchName() << "' was dropped in the first input file but is present in '"
                  << fileName << "'.\n";
    }
    return differences.str();
  }
}  // namespace edm
