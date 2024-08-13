#pragma once

#include <chrono>
#include <string>
#include <string_view>
#include <type_traits>

namespace precice::profiling {

/// Tag to annotate fundamental events
struct FundamentalTag {
};

/// Tag to annotate synchronized events
struct SynchronizeTag {
};

/// Convenience instance of the @ref FundamentalTag
static constexpr FundamentalTag Fundamental{};
/// Convenience instance of the @ref SynchronizeTag
static constexpr SynchronizeTag Synchronize{};

// Is the type a valid options tag?
template <typename T>
static constexpr bool isOptionsTag = std::is_same_v<T, FundamentalTag> || std::is_same_v<T, SynchronizeTag>;

/** Represents an event that can be started and stopped.
 *
 * Also allows to attach data in a key-value format using @ref addData()
 *
 * The event keeps minimal state. Events are passed to the @ref EventRegistry.
 */
class Event {
public:
  enum class State : bool {
    STOPPED = false,
    RUNNING = true
  };

  /// Default clock type. All other chrono types are derived from it.
  using Clock = std::chrono::steady_clock;

  /// Name used to identify the timer. Events of the same name are accumulated to
  std::string name;

  struct Options {
    bool synchronized = false;
    bool fundamental  = false;

    void handle(FundamentalTag)
    {
      fundamental = true;
    }
    void handle(SynchronizeTag)
    {
      synchronized = true;
    }
  };

  template <typename... Args>
  constexpr Options optionsFromTags(Args... args)
  {
    static_assert((isOptionsTag<Args> && ...), "The Event only accepts tags as arguments.");
    Options options;
    (options.handle(args), ...);
    return options;
  }

  template <typename... Args>
  Event(std::string_view eventName, Args... args)
      : Event(eventName, optionsFromTags(args...))
  {
  }

  Event(std::string_view eventName, Options options);

  Event(Event &&) = default;
  Event &operator=(Event &&) = default;

  // Copies would lead to duplicate entries
  Event(const Event &other) = delete;
  Event &operator=(const Event &) = delete;

  /// Stops the event if it's running and report its times to the EventRegistry
  ~Event();

  /// Starts or restarts a stopped event.
  void start();

  /// Stops a running event.
  void stop();

  /// Adds named integer data, associated to an event.
  void addData(std::string_view key, int value);

private:
  int   _eid;
  int   _sid{-1};
  State _state = State::STOPPED;
  bool  _fundamental{false};
  bool  _synchronize{false};
};

/// Class that changes the prefix in its scope
class ScopedEventPrefix {
public:
  ScopedEventPrefix(std::string_view name);

  ~ScopedEventPrefix();

  void pop();

private:
  std::string previousName = "";
};

} // namespace precice::profiling
