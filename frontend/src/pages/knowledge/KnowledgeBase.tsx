import { useState, useMemo } from 'react'
import { useNavigate } from 'react-router-dom'
import { Card, Row, Col, Typography, Input, Space, Button, Tag, List, Avatar, Tabs, Badge, Progress } from 'antd'
import {
  BookOutlined,
  SearchOutlined,
  FireOutlined,
  ClockCircleOutlined,
  LikeOutlined,
  ReadOutlined,
  RocketOutlined,
  BulbOutlined,
  QuestionCircleOutlined,
  CodeOutlined,
} from '@ant-design/icons'
import { PageHeader, PageSkeleton, FadeIn, EnhancedEmptyState } from '@/components/common'
import type { KnowledgeArticle, KnowledgeCategory, Tutorial } from '@/types/analytics'
import { usePopularArticles, useArticlesByCategory, useTutorials, useKnowledgeSearch } from '@/hooks/queries'
import { DesignTokens } from '@/styles/design-tokens'
import './KnowledgeBase.css'

const { Title, Text, Paragraph } = Typography
const { Search } = Input
const { TabPane } = Tabs

const categoryIcons: Record<KnowledgeCategory, any> = {
  'getting-started': <RocketOutlined />,
  pipelines: <CodeOutlined />,
  'data-management': <BookOutlined />,
  analysis: <BulbOutlined />,
  'quality-control': <FireOutlined />,
  troubleshooting: <QuestionCircleOutlined />,
  'best-practices': <LikeOutlined />,
  'api-reference': <CodeOutlined />,
  faq: <QuestionCircleOutlined />,
}

const categoryNames: Record<KnowledgeCategory, string> = {
  'getting-started': '快速开始',
  pipelines: '流程管道',
  'data-management': '数据管理',
  analysis: '数据分析',
  'quality-control': '质量控制',
  troubleshooting: '故障排除',
  'best-practices': '最佳实践',
  'api-reference': 'API参考',
  faq: '常见问题',
}

export const KnowledgeBase: React.FC = () => {
  const navigate = useNavigate()
  const [searchQuery, setSearchQuery] = useState('')
  const [selectedCategory, setSelectedCategory] = useState<KnowledgeCategory | 'all'>('all')

  const categoryFilter = selectedCategory !== 'all' ? selectedCategory : undefined

  // Search overrides category-based article browsing when present
  const isSearching = !!searchQuery
  const searchResults = useKnowledgeSearch(searchQuery, categoryFilter)
  const popularList = usePopularArticles(20)
  const categoryList = useArticlesByCategory(selectedCategory)
  const popularSidebar = usePopularArticles(3)
  const tutorialsQuery = useTutorials(categoryFilter)

  const articles: KnowledgeArticle[] = useMemo(() => {
    if (isSearching) {
      return (searchResults.data as KnowledgeArticle[]) ?? []
    }
    if (selectedCategory === 'all') {
      return (popularList.data as KnowledgeArticle[]) ?? []
    }
    return (categoryList.data as KnowledgeArticle[]) ?? []
  }, [isSearching, selectedCategory, searchResults.data, popularList.data, categoryList.data])

  const popularArticles: KnowledgeArticle[] = (popularSidebar.data as KnowledgeArticle[]) ?? []
  const tutorials: Tutorial[] = (tutorialsQuery.data as Tutorial[]) ?? []

  const loading =
    (isSearching && searchResults.isLoading) ||
    (!isSearching && selectedCategory === 'all' && popularList.isLoading) ||
    (!isSearching && selectedCategory !== 'all' && categoryList.isLoading) ||
    popularSidebar.isLoading ||
    tutorialsQuery.isLoading

  const initialLoad = articles.length === 0 && popularArticles.length === 0 && loading

  const handleSearch = (value: string) => {
    setSearchQuery(value)
  }

  const getDifficultyColor = (difficulty: string): string => {
    const colors = {
      beginner: 'green',
      intermediate: 'orange',
      advanced: 'red',
    }
    return colors[difficulty as keyof typeof colors] || 'blue'
  }

  const categories = Object.keys(categoryNames) as KnowledgeCategory[]

  // Show skeleton on initial load
  if (initialLoad && loading) {
    return <PageSkeleton hasHeader rows={6} />
  }

  return (
    <div className="knowledge-base">
      <FadeIn direction="up" delay={0} duration={300}>
        <PageHeader title="知识库" subtitle="教程、文档和最佳实践" icon={<BookOutlined />} />
      </FadeIn>

      {/* Search Bar */}
      <FadeIn direction="up" delay={50} duration={300}>
        <Card className="search-card" style={{ marginBottom: 24 }}>
          <Search
            placeholder="搜索文档、教程或问题..."
            allowClear
            enterButton={<SearchOutlined />}
            size="large"
            onSearch={handleSearch}
            style={{ maxWidth: 600 }}
          />
        </Card>
      </FadeIn>

      <Row gutter={[24, 24]}>
        {/* Left Sidebar - Categories */}
        <Col xs={24} md={6}>
          <FadeIn direction="up" delay={100} duration={300}>
            <Card title="分类" className="category-card">
              <Space direction="vertical" style={{ width: '100%' }} size="small">
                <Button
                  type={selectedCategory === 'all' ? 'primary' : 'text'}
                  block
                  onClick={() => setSelectedCategory('all')}
                >
                  全部文章
                </Button>
                {categories.map((category) => (
                  <Button
                    key={category}
                    type={selectedCategory === category ? 'primary' : 'text'}
                    icon={categoryIcons[category]}
                    block
                    onClick={() => setSelectedCategory(category)}
                  >
                    {categoryNames[category]}
                  </Button>
                ))}
              </Space>
            </Card>
          </FadeIn>

          {/* Popular Articles */}
          <FadeIn direction="up" delay={150} duration={300}>
            <Card
              title={
                <Space>
                  <FireOutlined style={{ color: DesignTokens.colors.error.main }} />
                  热门文章
                </Space>
              }
              className="popular-card"
              style={{ marginTop: 16 }}
            >
              <List
                dataSource={popularArticles}
                renderItem={(item) => (
                  <List.Item className="popular-item" onClick={() => navigate(`/knowledge/articles/${item.id}`)}>
                    <List.Item.Meta
                      title={
                        <Text strong style={{ fontSize: 13 }}>
                          {item.title}
                        </Text>
                      }
                      description={
                        <Space size="small">
                          <Text type="secondary" style={{ fontSize: 11 }}>
                            <ReadOutlined /> {item.views}
                          </Text>
                        </Space>
                      }
                    />
                  </List.Item>
                )}
              />
            </Card>
          </FadeIn>
        </Col>

        {/* Main Content */}
        <Col xs={24} md={18}>
          <FadeIn direction="up" delay={200} duration={300}>
            <Tabs defaultActiveKey="articles">
              <TabPane
                tab={
                  <Space>
                    <BookOutlined />
                    文章
                    <Badge count={articles.length} showZero />
                  </Space>
                }
                key="articles"
              >
                {articles.length === 0 ? (
                  <EnhancedEmptyState
                    type="noData"
                    title="No articles found"
                    description="没有找到相关文章"
                    size="default"
                  />
                ) : (
                  <List
                    dataSource={articles}
                    renderItem={(article) => (
                      <Card
                        className="article-card hover-lift"
                        hoverable
                        onClick={() => navigate(`/knowledge/articles/${article.id}`)}
                        style={{ marginBottom: 16 }}
                      >
                        <List.Item.Meta
                          avatar={
                            <Avatar
                              icon={categoryIcons[article.category]}
                              style={{ backgroundColor: DesignTokens.colors.primary.main }}
                            />
                          }
                          title={
                            <Space>
                              <Text strong style={{ fontSize: 16 }}>
                                {article.title}
                              </Text>
                              <Tag color={getDifficultyColor(article.difficulty)}>
                                {article.difficulty === 'beginner'
                                  ? '初级'
                                  : article.difficulty === 'intermediate'
                                    ? '中级'
                                    : '高级'}
                              </Tag>
                            </Space>
                          }
                          description={
                            <Space direction="vertical" style={{ width: '100%' }} size="small">
                              <Paragraph type="secondary" ellipsis={{ rows: 2 }}>
                                {article.summary}
                              </Paragraph>
                              <Space wrap>
                                {article.tags.map((tag) => (
                                  <Tag key={tag}>{tag}</Tag>
                                ))}
                              </Space>
                              <Space split="|" size="small">
                                <Text type="secondary" style={{ fontSize: 12 }}>
                                  <ClockCircleOutlined /> {article.estimatedReadTime} 分钟
                                </Text>
                                <Text type="secondary" style={{ fontSize: 12 }}>
                                  <ReadOutlined /> {article.views} 阅读
                                </Text>
                                <Text type="secondary" style={{ fontSize: 12 }}>
                                  <LikeOutlined /> {article.helpful} 有帮助
                                </Text>
                              </Space>
                            </Space>
                          }
                        />
                      </Card>
                    )}
                  />
                )}
              </TabPane>

              <TabPane
                tab={
                  <Space>
                    <RocketOutlined />
                    教程
                    <Badge count={tutorials.length} showZero />
                  </Space>
                }
                key="tutorials"
              >
                <Row gutter={[16, 16]}>
                  {tutorials.map((tutorial) => (
                    <Col xs={24} md={12} key={tutorial.id}>
                      <Card
                        className="tutorial-card hover-lift"
                        hoverable
                        onClick={() => navigate(`/knowledge/tutorials/${tutorial.id}`)}
                      >
                        <Space direction="vertical" style={{ width: '100%' }} size="middle">
                          <div>
                            <Title level={5} style={{ marginBottom: 4 }}>
                              {tutorial.title}
                            </Title>
                            <Space>
                              <Tag color={getDifficultyColor(tutorial.difficulty)}>
                                {tutorial.difficulty === 'beginner'
                                  ? '初级'
                                  : tutorial.difficulty === 'intermediate'
                                    ? '中级'
                                    : '高级'}
                              </Tag>
                              <Text type="secondary" style={{ fontSize: 12 }}>
                                <ClockCircleOutlined /> {tutorial.estimatedDuration} 分钟
                              </Text>
                            </Space>
                          </div>
                          <Paragraph type="secondary">{tutorial.description}</Paragraph>
                          <div>
                            <Text type="secondary" style={{ fontSize: 12 }}>
                              学习成果:
                            </Text>
                            <ul style={{ margin: '8px 0 0 20px', padding: 0 }}>
                              {tutorial.outcomes.map((outcome, idx) => (
                                <li key={idx}>
                                  <Text style={{ fontSize: 12 }}>{outcome}</Text>
                                </li>
                              ))}
                            </ul>
                          </div>
                          {tutorial.progress > 0 && (
                            <div>
                              <Text type="secondary" style={{ fontSize: 12 }}>
                                进度: {tutorial.progress}%
                              </Text>
                              <Progress percent={tutorial.progress} showInfo={false} />
                            </div>
                          )}
                        </Space>
                      </Card>
                    </Col>
                  ))}
                </Row>
              </TabPane>
            </Tabs>
          </FadeIn>
        </Col>
      </Row>
    </div>
  )
}
