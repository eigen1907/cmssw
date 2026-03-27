#ifndef RecoLocalMuon_RPCClusterPhase2_h
#define RecoLocalMuon_RPCClusterPhase2_h

class RPCClusterPhase2 {
public:
  RPCClusterPhase2();
  RPCClusterPhase2(int firstStrip, int lastStrip, int bx);
  ~RPCClusterPhase2();

  int firstStrip() const;
  int lastStrip() const;
  int clusterSize() const;
  int bx() const;

  bool hasTime() const;
  float time() const;
  float timeRMS() const;

  bool hasY() const;
  float y() const;
  float yRMS() const;

  void addTime(float time);
  void addY(float y);
  void merge(const RPCClusterPhase2& other);

  bool operator<(const RPCClusterPhase2& other) const;
  bool operator==(const RPCClusterPhase2& other) const;
  
  bool isAdjacent(const RPCClusterPhase2& other) const;
  bool isMaskedAdjacent(const RPCClusterPhase2& other) const;

private:
  int fstrip_;
  int lstrip_;
  int bx_;

  float sumTime_;
  float sumTime2_;
  unsigned int nTime_;

  float sumY_;
  float sumY2_;
  unsigned int nY_;
};

#endif