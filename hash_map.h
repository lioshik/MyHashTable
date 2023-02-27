// Короче это кукушка на двух массивах (пж не путать с кукушкой на одном массиве - это принципиально идейно два разных алгоритма не похожих друг на друга)
#pragma once

#include <cassert>
#include <list>
#include <memory>
#include <stdexcept>
#include <vector>

const size_t INITIAL_CAPACITY_1 = 40;
const size_t INITIAL_CAPACITY_2 = 41;
const double LOAD_FACTOR = 0.33;

template<typename KeyType, typename ValueType, typename Hasher = std::hash<KeyType>>
class HashMap {
private:
    std::vector<std::unique_ptr<std::pair<const KeyType, ValueType>>> array_1_, array_2_;
    std::list<std::unique_ptr<std::pair<const KeyType, ValueType>>> bad_items_;
    size_t item_count_;
    Hasher hasher_;

    void RebuildIfNeeded();

public:
    HashMap();

    HashMap(Hasher hasher);

    template<typename IteratorType>
    HashMap(IteratorType begin, IteratorType end);

    template<typename IteratorType>
    HashMap(IteratorType begin, IteratorType end, Hasher hasher);

    HashMap(std::initializer_list<std::pair<KeyType, ValueType>> given_items);

    HashMap(std::initializer_list<std::pair<KeyType, ValueType>> given_items, Hasher hasher);

    ~HashMap() = default;

    HashMap(const HashMap<KeyType, ValueType, Hasher> &other);

    HashMap &operator=(const HashMap<KeyType, ValueType, Hasher> &other);

    size_t size() const;

    size_t bad_size() const;

    bool empty() const;

    Hasher hash_function() const;

    void insert(std::pair<KeyType, ValueType> new_pair);

    void erase(KeyType key);

    class iterator {
    private:
        friend class HashMap<KeyType, ValueType, Hasher>;

        enum IteratingState {
            FIRST_ARRAY,
            SECOND_ARRAY,
            BAD_ITEMS
        };

        HashMap<KeyType, ValueType, Hasher> *original_map_;
        IteratingState iteration_stage_;
        size_t index_;
        typename std::list<std::unique_ptr<std::pair<const KeyType, ValueType>>>::iterator bad_items_iterator_;

        void move_to_nearest_item();

        void move_to_next_item();

    public:
        iterator();

        iterator(HashMap<KeyType, ValueType, Hasher> *original_map, IteratingState iteration_stage, size_t index,
                 typename std::list<std::unique_ptr<std::pair<const KeyType, ValueType>>>::iterator bad_items_iterator);

        std::pair<const KeyType, ValueType> &operator*() const;

        std::pair<const KeyType, ValueType> *operator->() const;

        iterator &operator++();

        iterator operator++(int);

        bool operator==(const iterator &other) const;

        bool operator!=(const iterator &other) const;
    };

    class const_iterator {
    private:
        friend class HashMap<KeyType, ValueType, Hasher>;

        enum IteratingState {
            FIRST_ARRAY,
            SECOND_ARRAY,
            BAD_ITEMS
        };

        const HashMap<KeyType, ValueType, Hasher> *original_map_;
        IteratingState iteration_stage_;
        size_t index_;
        typename std::list<std::unique_ptr<std::pair<const KeyType, ValueType>>>::const_iterator bad_items_iterator_;

        void move_to_nearest_item();

        void move_to_next_item();

    public:
        const_iterator();

        const_iterator(const HashMap<KeyType, ValueType, Hasher> *original_map, IteratingState iteration_stage,
                       size_t index,
                       typename std::list<std::unique_ptr<std::pair<const KeyType, ValueType>>>::const_iterator bad_items_iterator);

        const std::pair<const KeyType, ValueType> &operator*() const;

        const std::pair<const KeyType, ValueType> *operator->() const;

        const_iterator &operator++();

        const_iterator operator++(int);

        bool operator==(const const_iterator &other) const;

        bool operator!=(const const_iterator &other) const;
    };

    iterator begin();

    const_iterator begin() const;

    iterator end();

    const_iterator end() const;

    iterator find(KeyType key);

    const_iterator find(KeyType key) const;

    ValueType &operator[](KeyType key);

    const ValueType &at(KeyType key) const;

    void clear();
};

///////////// HashMap

template<typename KeyType, typename ValueType, typename Hasher>
void HashMap<KeyType, ValueType, Hasher>::RebuildIfNeeded() {
  if ((double) (array_1_.size() + array_2_.size()) * LOAD_FACTOR > item_count_) {
    return;
  }
  std::vector<std::unique_ptr<std::pair<const KeyType, ValueType>>> buffer;
  buffer.reserve(item_count_);

  for (auto &item: array_1_) {
    if (item.get() != nullptr) {
      buffer.push_back(std::move(item));
    }
  }
  for (auto &item: array_2_) {
    if (item.get() != nullptr) {
      buffer.push_back(std::move(item));
    }
  }
  for (auto &item: bad_items_) {
    if (item.get() != nullptr) {
      buffer.push_back(std::move(item));
    }
  }

  array_1_.resize(array_1_.size() * 2 + 1);
  array_2_.resize(array_2_.size() * 2 + 1);
  bad_items_.clear();
  for (auto &item: array_1_) {
    item = nullptr;
  }
  for (auto &item: array_2_) {
    item = nullptr;
  }
  item_count_ = 0;

  for (const auto &item: buffer) {
    insert(std::pair<KeyType, ValueType>(item->first, item->second));
  }
}

template<typename KeyType, typename ValueType, typename Hasher>
HashMap<KeyType, ValueType, Hasher>::HashMap() : array_1_(INITIAL_CAPACITY_1), array_2_(INITIAL_CAPACITY_2),
                                                 bad_items_(), item_count_(0), hasher_() {
}

template<typename KeyType, typename ValueType, typename Hasher>
HashMap<KeyType, ValueType, Hasher>::HashMap(Hasher hasher) : array_1_(INITIAL_CAPACITY_1),
                                                              array_2_(INITIAL_CAPACITY_2),
                                                              bad_items_(), item_count_(0), hasher_(hasher) {
}

template<typename KeyType, typename ValueType, typename Hasher>
template<typename IteratorType>
HashMap<KeyType, ValueType, Hasher>::HashMap(IteratorType begin, IteratorType end)
        : array_1_(INITIAL_CAPACITY_1), array_2_(INITIAL_CAPACITY_2), bad_items_(), item_count_(0), hasher_() {
  for (IteratorType iterator = begin; iterator != end; ++iterator) {
    insert(std::pair<KeyType, ValueType>(iterator->first, iterator->second));
  }
}

template<typename KeyType, typename ValueType, typename Hasher>
template<typename IteratorType>
HashMap<KeyType, ValueType, Hasher>::HashMap(IteratorType begin, IteratorType end, Hasher hasher)
        : array_1_(INITIAL_CAPACITY_1), array_2_(INITIAL_CAPACITY_2), bad_items_(), item_count_(0), hasher_(hasher) {
  for (IteratorType iterator = begin; iterator != end; ++iterator) {
    insert(std::pair<KeyType, ValueType>(iterator->first, iterator->second));
  }
}

template<typename KeyType, typename ValueType, typename Hasher>
HashMap<KeyType, ValueType, Hasher>::HashMap(std::initializer_list<std::pair<KeyType, ValueType>> given_items)
        : array_1_(INITIAL_CAPACITY_1), array_2_(INITIAL_CAPACITY_2), bad_items_(), item_count_(0), hasher_() {
  for (const auto &item: given_items) {
    insert(std::pair<KeyType, ValueType>(item.first, item.second));
  }
}

template<typename KeyType, typename ValueType, typename Hasher>
HashMap<KeyType, ValueType, Hasher>::HashMap(std::initializer_list<std::pair<KeyType, ValueType>> given_items,
                                             Hasher hasher) : array_1_(INITIAL_CAPACITY_1),
                                                              array_2_(INITIAL_CAPACITY_2), bad_items_(),
                                                              item_count_(0), hasher_(hasher) {
  for (const auto &item: given_items) {
    insert(std::pair<KeyType, ValueType>(item.first, item.second));
  }
}

template<typename KeyType, typename ValueType, typename Hasher>
size_t HashMap<KeyType, ValueType, Hasher>::size() const {
  return item_count_;
}

template<typename KeyType, typename ValueType, typename Hasher>
size_t HashMap<KeyType, ValueType, Hasher>::bad_size() const {
  return bad_items_.size();
}

template<typename KeyType, typename ValueType, typename Hasher>
bool HashMap<KeyType, ValueType, Hasher>::empty() const {
  return size() == 0;
}

template<typename KeyType, typename ValueType, typename Hasher>
Hasher HashMap<KeyType, ValueType, Hasher>::hash_function() const {
  return hasher_;
}

template<typename KeyType, typename ValueType, typename Hasher>
void HashMap<KeyType, ValueType, Hasher>::insert(std::pair<KeyType, ValueType> new_pair) {
  if (find(new_pair.first) != end()) {
    return;
  }
  RebuildIfNeeded();

  size_t initial_hash_1 = hasher_(new_pair.first) % array_1_.size();
  size_t initial_hash_2 = hasher_(new_pair.first) % array_2_.size();

  if (array_1_[initial_hash_1] == nullptr) {
    ++item_count_;
    array_1_[initial_hash_1] = std::make_unique<std::pair<const KeyType, ValueType>>(new_pair.first, new_pair.second);
    return;
  }

  if (array_2_[initial_hash_2] == nullptr) {
    ++item_count_;
    array_2_[initial_hash_2] = std::make_unique<std::pair<const KeyType, ValueType>>(new_pair.first, new_pair.second);
    return;
  }

  auto buffer = std::make_unique<std::pair<const KeyType, ValueType>>(new_pair.first, new_pair.second);
  auto current_hash_1 = initial_hash_1;
  auto current_hash_2 = hasher_(array_1_[initial_hash_1]->first) % array_2_.size();
  std::swap(buffer, array_1_[initial_hash_1]);

  size_t iteration_count = 0;
  while (!(current_hash_1 == initial_hash_1 && current_hash_2 == initial_hash_2)) {
    if (array_1_[current_hash_1] == nullptr) {
      std::swap(array_1_[current_hash_1], buffer);
      ++item_count_;
      return;
    }
    if (array_2_[current_hash_2] == nullptr) {
      std::swap(array_2_[current_hash_2], buffer);
      ++item_count_;
      return;
    }

    if (iteration_count % 2 == 0) {
      current_hash_1 = hasher_(array_2_[current_hash_2]->first) % array_1_.size();
      std::swap(buffer, array_2_[current_hash_2]);
    } else {
      current_hash_2 = hasher_(array_1_[current_hash_1]->first) % array_2_.size();
      std::swap(buffer, array_1_[current_hash_1]);
    }
    ++iteration_count;
  }

  if (current_hash_1 == initial_hash_1 && current_hash_2 == initial_hash_2) {
    bad_items_.push_front(std::move(buffer));
    ++item_count_;
    return;
  }
  assert(false); // we shouldn't be here
}

template<typename KeyType, typename ValueType, typename Hasher>
void HashMap<KeyType, ValueType, Hasher>::erase(KeyType key) {
  size_t hash1 = hasher_(key) % array_1_.size();
  size_t hash2 = hasher_(key) % array_2_.size();
  if (array_1_[hash1] && array_1_[hash1]->first == key) {
    array_1_[hash1] = nullptr;
    --item_count_;
    return;
  }

  if (array_2_[hash2] && array_2_[hash2]->first == key) {
    array_2_[hash2] = nullptr;
    --item_count_;
    return;
  }

  for (auto item_iterator = bad_items_.begin(); item_iterator != bad_items_.end(); ++item_iterator) {
    if ((*item_iterator)->first == key) {
      *item_iterator = nullptr;
      bad_items_.erase(item_iterator);
      --item_count_;
      return;
    }
  }
}

template<typename KeyType, typename ValueType, typename Hasher>
typename HashMap<KeyType, ValueType, Hasher>::iterator HashMap<KeyType, ValueType, Hasher>::begin() {
  return HashMap::iterator(this, iterator::IteratingState::FIRST_ARRAY, 0, bad_items_.begin());
}

template<typename KeyType, typename ValueType, typename Hasher>
typename HashMap<KeyType, ValueType, Hasher>::const_iterator HashMap<KeyType, ValueType, Hasher>::begin() const {
  return HashMap::const_iterator(this, const_iterator::IteratingState::FIRST_ARRAY, 0, bad_items_.begin());
}

template<typename KeyType, typename ValueType, typename Hasher>
typename HashMap<KeyType, ValueType, Hasher>::iterator HashMap<KeyType, ValueType, Hasher>::end() {
  return HashMap::iterator(this, iterator::IteratingState::BAD_ITEMS, 0, bad_items_.end());
}

template<typename KeyType, typename ValueType, typename Hasher>
typename HashMap<KeyType, ValueType, Hasher>::const_iterator HashMap<KeyType, ValueType, Hasher>::end() const {
  return HashMap::const_iterator(this, const_iterator::IteratingState::BAD_ITEMS, 0, bad_items_.end());
}

template<typename KeyType, typename ValueType, typename Hasher>
typename HashMap<KeyType, ValueType, Hasher>::iterator HashMap<KeyType, ValueType, Hasher>::find(KeyType key) {
  size_t hash1 = hasher_(key) % array_1_.size();
  size_t hash2 = hasher_(key) % array_2_.size();
  if (array_1_[hash1] && array_1_[hash1]->first == key) {
    return iterator(this, iterator::IteratingState::FIRST_ARRAY, hash1, bad_items_.begin());
  }

  if (array_2_[hash2] && array_2_[hash2]->first == key) {
    return iterator(this, iterator::IteratingState::SECOND_ARRAY, hash2, bad_items_.begin());
  }

  for (auto item_iterator = bad_items_.begin(); item_iterator != bad_items_.end(); ++item_iterator) {
    if ((*item_iterator)->first == key) {
      return iterator(this, iterator::IteratingState::BAD_ITEMS, 0, item_iterator);
    }
  }

  return end();
}

template<typename KeyType, typename ValueType, typename Hasher>
typename HashMap<KeyType, ValueType, Hasher>::const_iterator
HashMap<KeyType, ValueType, Hasher>::find(KeyType key) const {
  size_t hash1 = hasher_(key) % array_1_.size();
  size_t hash2 = hasher_(key) % array_2_.size();
  if (array_1_[hash1] && array_1_[hash1]->first == key) {
    return const_iterator(this, const_iterator::IteratingState::FIRST_ARRAY, hash1, bad_items_.begin());
  }

  if (array_2_[hash2] && array_2_[hash2]->first == key) {
    return const_iterator(this, const_iterator::IteratingState::SECOND_ARRAY, hash2, bad_items_.begin());
  }

  for (auto item_iterator = bad_items_.begin(); item_iterator != bad_items_.end(); ++item_iterator) {
    if ((*item_iterator)->first == key) {
      return const_iterator(this, const_iterator::IteratingState::BAD_ITEMS, 0, item_iterator);
    }
  }

  return end();
}

template<typename KeyType, typename ValueType, typename Hasher>
ValueType &HashMap<KeyType, ValueType, Hasher>::operator[](KeyType key) {
  size_t hash1 = hasher_(key) % array_1_.size();
  size_t hash2 = hasher_(key) % array_2_.size();
  if (array_1_[hash1] && array_1_[hash1]->first == key) {
    return array_1_[hash1]->second;
  }

  if (array_2_[hash2] && array_2_[hash2]->first == key) {
    return array_2_[hash2]->second;
  }

  for (auto item_iterator = bad_items_.begin(); item_iterator != bad_items_.end(); ++item_iterator) {
    if ((*item_iterator)->first == key) {
      return (*item_iterator)->second;
    }
  }

  insert(std::pair<KeyType, ValueType>(key, ValueType()));
  return operator[](key);
}

template<typename KeyType, typename ValueType, typename Hasher>
const ValueType &HashMap<KeyType, ValueType, Hasher>::at(KeyType key) const {
  size_t hash1 = hasher_(key) % array_1_.size();
  size_t hash2 = hasher_(key) % array_2_.size();
  if (array_1_[hash1] && array_1_[hash1]->first == key) {
    return array_1_[hash1]->second;
  }

  if (array_2_[hash2] && array_2_[hash2]->first == key) {
    return array_2_[hash2]->second;
  }

  for (auto item_iterator = bad_items_.begin(); item_iterator != bad_items_.end(); ++item_iterator) {
    if ((*item_iterator)->first == key) {
      return (*item_iterator)->second;
    }
  }

  throw std::out_of_range("error: map does not contain key");
}

template<typename KeyType, typename ValueType, typename Hasher>
void HashMap<KeyType, ValueType, Hasher>::clear() {
  array_1_.clear();
  array_2_.clear();
  array_1_.resize(INITIAL_CAPACITY_1);
  array_2_.resize(INITIAL_CAPACITY_2);
  bad_items_.clear();
  item_count_ = 0;
}

template<typename KeyType, typename ValueType, typename Hasher>
HashMap<KeyType, ValueType, Hasher>::HashMap(const HashMap<KeyType, ValueType, Hasher> &other)
        : array_1_(other.array_1_.size()), array_2_(other.array_2_.size()), bad_items_(), item_count_(0),
          hasher_(other.hasher_) {
  for (const auto &item: other) {
    insert(std::pair<KeyType, ValueType>(item.first, item.second));
  }
}

template<typename KeyType, typename ValueType, typename Hasher>
HashMap<KeyType, ValueType, Hasher> &
HashMap<KeyType, ValueType, Hasher>::operator=(const HashMap<KeyType, ValueType, Hasher> &other) {
  if (this == &other) {
    return *this;
  }
  clear();
  array_1_.resize(other.array_1_.size());
  array_2_.resize(other.array_2_.size());
  item_count_ = other.item_count_;
  hasher_ = other.hasher_;
  for (size_t i = 0; i < array_1_.size(); ++i) {
    if (other.array_1_[i] != nullptr) {
      array_1_[i] = std::make_unique<std::pair<const KeyType, ValueType>>(other.array_1_[i]->first,
                                                                          other.array_1_[i]->second);
    }
  }
  for (size_t i = 0; i < array_2_.size(); ++i) {
    if (other.array_2_[i] != nullptr) {
      array_2_[i] = std::make_unique<std::pair<const KeyType, ValueType>>(other.array_2_[i]->first,
                                                                          other.array_2_[i]->second);
    }
  }
  for (const auto &item : other.bad_items_) {
    bad_items_.push_back(std::make_unique<std::pair<const KeyType, ValueType>>(item->first, item->second));
  }

  return *this;
}

/////////////////////////////////////////// iterator

template<typename KeyType, typename ValueType, typename Hasher>
HashMap<KeyType, ValueType, Hasher>::iterator::iterator() : original_map_(nullptr), iteration_stage_(FIRST_ARRAY),
                                                            index_(0), bad_items_iterator_() {
}

template<typename KeyType, typename ValueType, typename Hasher>
void HashMap<KeyType, ValueType, Hasher>::iterator::move_to_next_item() {
  switch (iteration_stage_) {
    case IteratingState::FIRST_ARRAY:
      ++index_;
      while (index_ < original_map_->array_1_.size() && original_map_->array_1_[index_] == nullptr) {
        ++index_;
      }
      if (index_ >= original_map_->array_1_.size() || original_map_->array_1_[index_] == nullptr) {
        index_ = 0;
        iteration_stage_ = IteratingState::SECOND_ARRAY;
        if (original_map_->array_2_[0] == nullptr) {
          move_to_nearest_item();
        }
      }
      break;
    case IteratingState::SECOND_ARRAY:
      ++index_;
      while (index_ < original_map_->array_2_.size() && original_map_->array_2_[index_] == nullptr) {
        ++index_;
      }
      if (index_ >= original_map_->array_2_.size() || original_map_->array_2_[index_] == nullptr) {
        index_ = 0;
        iteration_stage_ = IteratingState::BAD_ITEMS;
        bad_items_iterator_ = original_map_->bad_items_.begin();
      }
      break;
    case IteratingState::BAD_ITEMS:
      ++bad_items_iterator_;
      break;
    default:
      assert(false); // we shouldn't be here
  }
}

template<typename KeyType, typename ValueType, typename Hasher>
void HashMap<KeyType, ValueType, Hasher>::iterator::move_to_nearest_item() {
  switch (iteration_stage_) {
    case IteratingState::FIRST_ARRAY:
      if (original_map_->array_1_[index_] == nullptr) {
        move_to_next_item();
      }
      break;
    case IteratingState::SECOND_ARRAY:
      if (original_map_->array_2_[index_] == nullptr) {
        move_to_next_item();
      }
      break;
    case IteratingState::BAD_ITEMS:
      break;
    default:
      assert(false); // we shouldn't be here
  }
}

template<typename KeyType, typename ValueType, typename Hasher>
HashMap<KeyType, ValueType, Hasher>::iterator::iterator(
        HashMap<KeyType, ValueType, Hasher> *original_map,
        IteratingState iteration_stage, size_t index,
        typename std::list<std::unique_ptr<std::pair<const KeyType, ValueType>>>::iterator bad_items_iterator)
        : original_map_(original_map), iteration_stage_(iteration_stage), index_(index),
          bad_items_iterator_(bad_items_iterator) {
  move_to_nearest_item();
}

template<typename KeyType, typename ValueType, typename Hasher>
std::pair<const KeyType, ValueType> &HashMap<KeyType, ValueType, Hasher>::iterator::operator*() const {
  switch (iteration_stage_) {
    case IteratingState::FIRST_ARRAY:
      return *(original_map_->array_1_[index_]);
    case IteratingState::SECOND_ARRAY:
      return *(original_map_->array_2_[index_]);
    case IteratingState::BAD_ITEMS:
      return bad_items_iterator_->operator*();
    default:
      assert(false); // we shouldn't be here
  }
}

template<typename KeyType, typename ValueType, typename Hasher>
std::pair<const KeyType, ValueType> *HashMap<KeyType, ValueType, Hasher>::iterator::operator->() const {
  switch (iteration_stage_) {
    case IteratingState::FIRST_ARRAY:
      return original_map_->array_1_[index_].get();
    case IteratingState::SECOND_ARRAY:
      return original_map_->array_2_[index_].get();
    case IteratingState::BAD_ITEMS:
      return bad_items_iterator_->operator->();
    default:
      assert(false); // we shouldn't be here
  }
}

template<typename KeyType, typename ValueType, typename Hasher>
typename HashMap<KeyType, ValueType, Hasher>::iterator &HashMap<KeyType, ValueType, Hasher>::iterator::operator++() {
  move_to_next_item();
  return *this;
}

template<typename KeyType, typename ValueType, typename Hasher>
typename HashMap<KeyType, ValueType, Hasher>::iterator HashMap<KeyType, ValueType, Hasher>::iterator::operator++(int) {
  auto result = iterator(*this);
  move_to_next_item();
  return result;
}

template<typename KeyType, typename ValueType, typename Hasher>
bool HashMap<KeyType, ValueType, Hasher>::iterator::operator==(const HashMap::iterator &other) const {
  assert(original_map_ != nullptr);
  assert(other.original_map_ != nullptr);
  if (original_map_ != other.original_map_) {
    return false;
  }

  if (iteration_stage_ != other.iteration_stage_) {
    return false;
  }
  switch (iteration_stage_) {
    case IteratingState::FIRST_ARRAY:
      return index_ == other.index_;
    case IteratingState::SECOND_ARRAY:
      return index_ == other.index_;
    case IteratingState::BAD_ITEMS:
      return bad_items_iterator_ == other.bad_items_iterator_;
    default:
      assert(false); // we shouldn't be here
  }
}

template<typename KeyType, typename ValueType, typename Hasher>
bool HashMap<KeyType, ValueType, Hasher>::iterator::operator!=(const HashMap::iterator &other) const {
  return !(this->operator==(other));
}

/////////////////////////////////////////// const_iterator

template<typename KeyType, typename ValueType, typename Hasher>
void HashMap<KeyType, ValueType, Hasher>::const_iterator::move_to_nearest_item() {
  switch (iteration_stage_) {
    case IteratingState::FIRST_ARRAY:
      if (original_map_->array_1_[index_] == nullptr) {
        move_to_next_item();
      }
      break;
    case IteratingState::SECOND_ARRAY:
      if (original_map_->array_2_[index_] == nullptr) {
        move_to_next_item();
      }
      break;
    case IteratingState::BAD_ITEMS:
      break;
    default:
      assert(false); // we shouldn't be here
  }
}

template<typename KeyType, typename ValueType, typename Hasher>
void HashMap<KeyType, ValueType, Hasher>::const_iterator::move_to_next_item() {
  switch (iteration_stage_) {
    case IteratingState::FIRST_ARRAY:
      ++index_;
      while (index_ < original_map_->array_1_.size() && original_map_->array_1_[index_] == nullptr) {
        ++index_;
      }
      if (index_ >= original_map_->array_1_.size() || original_map_->array_1_[index_] == nullptr) {
        index_ = 0;
        iteration_stage_ = IteratingState::SECOND_ARRAY;
        if (original_map_->array_2_[0] == nullptr) {
          move_to_nearest_item();
        }
      }
      break;
    case IteratingState::SECOND_ARRAY:
      ++index_;
      while (index_ < original_map_->array_2_.size() && original_map_->array_2_[index_] == nullptr) {
        ++index_;
      }
      if (index_ >= original_map_->array_2_.size() || original_map_->array_2_[index_] == nullptr) {
        index_ = 0;
        iteration_stage_ = IteratingState::BAD_ITEMS;
        bad_items_iterator_ = original_map_->bad_items_.begin();
      }
      break;
    case IteratingState::BAD_ITEMS:
      ++bad_items_iterator_;
      break;
    default:
      assert(false); // we shouldn't be here
  }
}

template<typename KeyType, typename ValueType, typename Hasher>
HashMap<KeyType, ValueType, Hasher>::const_iterator::const_iterator()
        : original_map_(nullptr), iteration_stage_(FIRST_ARRAY), index_(0), bad_items_iterator_() {
}

template<typename KeyType, typename ValueType, typename Hasher>
HashMap<KeyType, ValueType, Hasher>::const_iterator::const_iterator(
        const HashMap<KeyType, ValueType, Hasher> *original_map, IteratingState
iteration_stage, size_t
        index,
        typename std::list<std::unique_ptr<std::pair<const KeyType, ValueType>>>::const_iterator
        bad_items_iterator)
        : original_map_(original_map), iteration_stage_(iteration_stage), index_(index),
          bad_items_iterator_(bad_items_iterator) {
  move_to_nearest_item();
}

template<typename KeyType, typename ValueType, typename Hasher>
const std::pair<const KeyType, ValueType> &HashMap<KeyType, ValueType, Hasher>::const_iterator::operator*() const {
  switch (iteration_stage_) {
    case IteratingState::FIRST_ARRAY:
      return *(original_map_->array_1_[index_]);
    case IteratingState::SECOND_ARRAY:
      return *(original_map_->array_2_[index_]);
    case IteratingState::BAD_ITEMS:
      return bad_items_iterator_->operator*();
    default:
      assert(false); // we shouldn't be here
  }
}

template<typename KeyType, typename ValueType, typename Hasher>
const std::pair<const KeyType, ValueType> *HashMap<KeyType, ValueType, Hasher>::const_iterator::operator->() const {
  switch (iteration_stage_) {
    case IteratingState::FIRST_ARRAY:
      return original_map_->array_1_[index_].get();
    case IteratingState::SECOND_ARRAY:
      return original_map_->array_2_[index_].get();
    case IteratingState::BAD_ITEMS:
      return bad_items_iterator_->operator->();
    default:
      assert(false); // we shouldn't be here
  }
}

template<typename KeyType, typename ValueType, typename Hasher>
typename HashMap<KeyType, ValueType, Hasher>::const_iterator &
HashMap<KeyType, ValueType, Hasher>::const_iterator::operator++() {
  move_to_next_item();
  return *this;
}

template<typename KeyType, typename ValueType, typename Hasher>
typename HashMap<KeyType, ValueType, Hasher>::const_iterator
HashMap<KeyType, ValueType, Hasher>::const_iterator::operator++(int) {
  auto result = const_iterator(*this);
  move_to_next_item();
  return result;
}

template<typename KeyType, typename ValueType, typename Hasher>
bool HashMap<KeyType, ValueType, Hasher>::const_iterator::operator==(const HashMap::const_iterator &other) const {
  assert(original_map_ != nullptr);
  assert(other.original_map_ != nullptr);
  if (original_map_ != other.original_map_) {
    return false;
  }

  if (iteration_stage_ != other.iteration_stage_) {
    return false;
  }
  switch (iteration_stage_) {
    case IteratingState::FIRST_ARRAY:
      return index_ == other.index_;
    case IteratingState::SECOND_ARRAY:
      return index_ == other.index_;
    case IteratingState::BAD_ITEMS:
      return bad_items_iterator_ == other.bad_items_iterator_;
    default:
      assert(false); // we shouldn't be here
  }
}

template<typename KeyType, typename ValueType, typename Hasher>
bool HashMap<KeyType, ValueType, Hasher>::const_iterator::operator!=(const HashMap::const_iterator &other) const {
  return !(this->operator==(other));
}

