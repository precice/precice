#ifndef GLOBALS_HPP_
#define GLOBALS_HPP_

/**
 * @brief For printing to the command line.
 *
 * @param message  An input stream such as: "Blabla = " << variable << ", ..."
 */
#ifdef STRUCTURE_DEBUG_MODE
#define STRUCTURE_DEBUG(message) \
   { \
      std::ostringstream conv; \
      conv << message; \
      std::cout << conv.str() << '\n'; \
   }
#else
#define STRUCTURE_DEBUG(message)
#endif

#define STRUCTURE_INFO(message) \
   { \
      std::ostringstream conv; \
      conv << message; \
      std::cout << conv.str() << '\n'; \
   }

#endif /* GLOBALS_HPP_ */
