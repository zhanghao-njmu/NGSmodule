import { useState, useEffect } from 'react'
import { useNavigate } from 'react-router-dom'
import {
  Card,
  Row,
  Col,
  Typography,
  Input,
  Space,
  Button,
  Tag,
  List,
  Avatar,
  Empty,
  Tabs,
  Badge,
  Spin,
  Progress,
} from 'antd'
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
import { PageHeader } from '@/components/common'
import type { KnowledgeArticle, KnowledgeCategory, Tutorial } from '@/types/analytics'
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
  const [loading, setLoading] = useState(false)
  const [_searchQuery, setSearchQuery] = useState('')
  const [selectedCategory, setSelectedCategory] = useState<KnowledgeCategory | 'all'>('all')
  const [articles, setArticles] = useState<KnowledgeArticle[]>([])
  const [popularArticles, setPopularArticles] = useState<KnowledgeArticle[]>([])
  const [tutorials, setTutorials] = useState<Tutorial[]>([])

  useEffect(() => {
    fetchKnowledgeContent()
  }, [selectedCategory])

  const fetchKnowledgeContent = async () => {
    setLoading(true)
    try {
      // TODO: Replace with actual API calls
      // const articles = await analyticsService.getArticlesByCategory(selectedCategory)
      // const popular = await analyticsService.getPopularArticles()
      // const tutorials = await analyticsService.getTutorials()

      // Mock data
      const mockArticles: KnowledgeArticle[] = [
        {
          id: '1',
          title: 'NGSmodule 快速入门指南',
          category: 'getting-started',
          tags: ['新手', '教程', '基础'],
          content: '详细的快速入门指南...',
          summary: '了解如何快速开始使用NGSmodule平台进行生物信息学分析',
          author: 'Admin',
          createdAt: '2024-01-15',
          updatedAt: '2024-03-20',
          views: 1245,
          helpful: 98,
          notHelpful: 2,
          difficulty: 'beginner',
          estimatedReadTime: 10,
          relatedArticles: ['2', '3'],
        },
        {
          id: '2',
          title: '如何配置RNA-seq分析流程',
          category: 'pipelines',
          tags: ['RNA-seq', '流程', '配置'],
          content: 'RNA-seq流程配置详解...',
          summary: '学习如何配置和运行RNA-seq分析流程',
          author: 'Expert',
          createdAt: '2024-02-01',
          updatedAt: '2024-03-18',
          views: 892,
          helpful: 76,
          notHelpful: 4,
          difficulty: 'intermediate',
          estimatedReadTime: 15,
          relatedArticles: ['1', '4'],
        },
        {
          id: '3',
          title: '数据质量控制最佳实践',
          category: 'quality-control',
          tags: ['QC', '质量', '最佳实践'],
          content: 'QC最佳实践详解...',
          summary: '掌握数据质量控制的关键步骤和最佳实践',
          author: 'Expert',
          createdAt: '2024-02-10',
          updatedAt: '2024-03-15',
          views: 756,
          helpful: 65,
          notHelpful: 3,
          difficulty: 'intermediate',
          estimatedReadTime: 12,
          relatedArticles: ['2', '5'],
        },
        {
          id: '4',
          title: '常见问题解答',
          category: 'faq',
          tags: ['FAQ', '问题', '解答'],
          content: '常见问题集合...',
          summary: '查找常见问题的答案',
          author: 'Admin',
          createdAt: '2024-01-20',
          updatedAt: '2024-03-22',
          views: 2134,
          helpful: 156,
          notHelpful: 12,
          difficulty: 'beginner',
          estimatedReadTime: 8,
          relatedArticles: ['1'],
        },
      ]

      setArticles(mockArticles)
      setPopularArticles(mockArticles.slice(0, 3))

      const mockTutorials: Tutorial[] = [
        {
          id: '1',
          title: '从零开始的RNA-seq分析',
          description: '完整的RNA-seq分析教程，从数据上传到结果解读',
          category: 'getting-started',
          difficulty: 'beginner',
          estimatedDuration: 45,
          steps: [],
          outcomes: ['掌握RNA-seq基本流程', '理解质量控制', '学会结果解读'],
          completed: false,
          progress: 0,
        },
        {
          id: '2',
          title: 'AI功能使用指南',
          description: '学习如何使用平台的AI智能功能提高工作效率',
          category: 'best-practices',
          difficulty: 'intermediate',
          estimatedDuration: 30,
          steps: [],
          outcomes: ['掌握AI推荐系统', '使用自动QC', '智能分组功能'],
          completed: false,
          progress: 0,
        },
      ]

      setTutorials(mockTutorials)
    } catch (error) {
      console.error('Failed to fetch knowledge content:', error)
    } finally {
      setLoading(false)
    }
  }

  const handleSearch = async (value: string) => {
    setSearchQuery(value)
    if (value) {
      setLoading(true)
      try {
        // TODO: Call search API
        // const results = await analyticsService.searchKnowledge(value, selectedCategory !== 'all' ? selectedCategory : undefined)
        // setArticles(results)
      } catch (error) {
        console.error('Search failed:', error)
      } finally {
        setLoading(false)
      }
    } else {
      fetchKnowledgeContent()
    }
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

  return (
    <div className="knowledge-base">
      <PageHeader title="知识库" subtitle="教程、文档和最佳实践" icon={<BookOutlined />} />

      {/* Search Bar */}
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

      <Row gutter={[24, 24]}>
        {/* Left Sidebar - Categories */}
        <Col xs={24} md={6}>
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

          {/* Popular Articles */}
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
        </Col>

        {/* Main Content */}
        <Col xs={24} md={18}>
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
              {loading ? (
                <div style={{ textAlign: 'center', padding: '48px' }}>
                  <Spin size="large" />
                </div>
              ) : articles.length === 0 ? (
                <Empty description="没有找到相关文章" />
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
        </Col>
      </Row>
    </div>
  )
}
