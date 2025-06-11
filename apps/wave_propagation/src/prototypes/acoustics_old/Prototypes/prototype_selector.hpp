
//  Contributions by Omar Durán and Romain Mottier

#ifndef prototype_selector_hpp
#define prototype_selector_hpp

void print_prototype_description(){
}

int prototype_selector(char **argv, EAcousticPrototype prototype){
    
    switch (prototype) {

      case EAcousticPrototype::OneFieldConvTest: {
	EllipticOneFieldConvergenceTest(argv);
      }
      break;

      case EAcousticPrototype::TwoFieldsConvTest: {
	EllipticTwoFieldsConvergenceTest(argv);
      }
      break;

      case EAcousticPrototype::IOneFieldAcoustic: {
        IHHOSecondOrder(argv);
      }
      break;

      case EAcousticPrototype::ITwoFieldsAcoustic: {
        IHHOFirstOrder(argv);
      }
      break;

      case EAcousticPrototype::ETwoFieldsAcoustic: {
        EHHOFirstOrder(argv);
      }
      break;

      case EAcousticPrototype::Hete1DIOneFieldAcoustic: {
        HeterogeneousIHHOSecondOrder(argv);
      }
      break;

      case EAcousticPrototype::Hete1DITwoFieldsAcoustic: {
        HeterogeneousIHHOFirstOrder(argv);
      }
      break;

      case EAcousticPrototype::Hete1DETwoFieldsAcoustic: {
        HeterogeneousEHHOFirstOrder(argv);
      }
      break;

      case EAcousticPrototype::Hete2DIOneFieldAcoustic: {
        HeterogeneousPulseIHHOSecondOrder(argv);
      }
      break;

      case EAcousticPrototype::Hete2DITwoFieldsAcoustic: {
        HeterogeneousPulseIHHOFirstOrder(argv);
      }
      break;

      case EAcousticPrototype::Hete2DETwoFieldsAcoustic: {
        HeterogeneousPulseEHHOFirstOrder(argv);
      }
      break;

      default: {
        throw std::invalid_argument("Prototype not available.");
      }
      break;

    }
    
}

#endif 
