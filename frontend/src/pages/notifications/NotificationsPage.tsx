import { useState, useEffect } from 'react'
import { useNavigate } from 'react-router-dom'
import {
  Card,
  List,
  Typography,
  Button,
  Space,
  Tag,
  Empty,
  Spin,
  Select,
  Checkbox,
  Pagination,
  Modal,
  message,
  Row,
  Col,
} from 'antd'
import {
  BellOutlined,
  CheckCircleOutlined,
  InfoCircleOutlined,
  WarningOutlined,
  CloseCircleOutlined,
  CheckOutlined,
  DeleteOutlined,
  FilterOutlined,
  InboxOutlined,
} from '@ant-design/icons'
import { PageHeader, StatisticCard } from '@/components/common'
import type { Notification, NotificationType, NotificationCategory } from '@/types/notification'
import { DesignTokens } from '@/styles/design-tokens'
import './NotificationsPage.css'

const { Text, Paragraph } = Typography
const { Option } = Select

export const NotificationsPage: React.FC = () => {
  const navigate = useNavigate()
  const [notifications, setNotifications] = useState<Notification[]>([])
  const [loading, setLoading] = useState(false)
  const [selectedIds, setSelectedIds] = useState<string[]>([])
  const [currentPage, setCurrentPage] = useState(1)
  const [pageSize, setPageSize] = useState(10)
  const [total, setTotal] = useState(0)
  const [filterType, setFilterType] = useState<NotificationType | 'all'>('all')
  const [filterCategory, setFilterCategory] = useState<NotificationCategory | 'all'>('all')
  const [filterRead, setFilterRead] = useState<'all' | 'read' | 'unread'>('all')

  // Statistics
  const [stats, setStats] = useState({
    total: 0,
    unread: 0,
    success: 0,
    error: 0,
    warning: 0,
    info: 0,
  })

  useEffect(() => {
    fetchNotifications()
    fetchStats()
  }, [currentPage, pageSize, filterType, filterCategory, filterRead])

  const fetchNotifications = async () => {
    setLoading(true)
    try {
      // TODO: Replace with actual API call
      // const response = await notificationService.getNotifications({
      //   page: currentPage,
      //   pageSize,
      //   filter: {
      //     type: filterType !== 'all' ? [filterType] : undefined,
      //     category: filterCategory !== 'all' ? [filterCategory] : undefined,
      //     read: filterRead === 'all' ? undefined : filterRead === 'read',
      //   },
      // })

      // Mock data
      const mockNotifications: Notification[] = [
        {
          id: '1',
          type: 'success',
          category: 'pipeline',
          title: 'Pipeline completed successfully',
          message: 'RNA-seq analysis pipeline has finished processing sample batch A',
          timestamp: new Date(Date.now() - 2 * 60 * 60 * 1000).toISOString(),
          read: false,
          actionUrl: '/results/123',
          actionText: 'View Results',
        },
        {
          id: '2',
          type: 'error',
          category: 'task',
          title: 'Task failed',
          message: 'Alignment task encountered an error during execution',
          timestamp: new Date(Date.now() - 5 * 60 * 60 * 1000).toISOString(),
          read: false,
          actionUrl: '/tasks/456',
          actionText: 'View Details',
        },
        {
          id: '3',
          type: 'info',
          category: 'system',
          title: 'System maintenance scheduled',
          message: 'Scheduled maintenance on Sunday 2:00 AM - 4:00 AM',
          timestamp: new Date(Date.now() - 24 * 60 * 60 * 1000).toISOString(),
          read: true,
        },
        {
          id: '4',
          type: 'warning',
          category: 'project',
          title: 'Storage limit warning',
          message: 'Project "Cancer Genomics" storage is 85% full',
          timestamp: new Date(Date.now() - 2 * 24 * 60 * 60 * 1000).toISOString(),
          read: true,
          actionUrl: '/items/789',
          actionText: 'Manage Storage',
        },
        {
          id: '5',
          type: 'success',
          category: 'sample',
          title: 'Samples uploaded',
          message: '24 new samples have been successfully uploaded',
          timestamp: new Date(Date.now() - 3 * 24 * 60 * 60 * 1000).toISOString(),
          read: true,
          actionUrl: '/samples',
          actionText: 'View Samples',
        },
      ]

      setNotifications(mockNotifications)
      setTotal(mockNotifications.length)
    } catch (error) {
      message.error('Failed to load notifications')
    } finally {
      setLoading(false)
    }
  }

  const fetchStats = async () => {
    try {
      // TODO: Replace with actual API call
      // const statsData = await notificationService.getStats()
      setStats({
        total: 45,
        unread: 12,
        success: 20,
        error: 5,
        warning: 8,
        info: 12,
      })
    } catch (error) {
      console.error('Failed to load stats:', error)
    }
  }

  const handleSelectAll = (checked: boolean) => {
    if (checked) {
      setSelectedIds(notifications.map((n) => n.id))
    } else {
      setSelectedIds([])
    }
  }

  const handleSelect = (id: string, checked: boolean) => {
    if (checked) {
      setSelectedIds([...selectedIds, id])
    } else {
      setSelectedIds(selectedIds.filter((sid) => sid !== id))
    }
  }

  const handleMarkAsRead = async (ids: string[]) => {
    try {
      // TODO: Call API
      // await notificationService.markMultipleAsRead(ids)
      setNotifications(
        notifications.map((n) => (ids.includes(n.id) ? { ...n, read: true } : n))
      )
      setSelectedIds([])
      message.success(`${ids.length} notification(s) marked as read`)
      fetchStats()
    } catch (error) {
      message.error('Failed to mark notifications as read')
    }
  }

  const handleMarkAllAsRead = async () => {
    try {
      // TODO: Call API
      // await notificationService.markAllAsRead()
      setNotifications(notifications.map((n) => ({ ...n, read: true })))
      message.success('All notifications marked as read')
      fetchStats()
    } catch (error) {
      message.error('Failed to mark all as read')
    }
  }

  const handleDelete = async (ids: string[]) => {
    Modal.confirm({
      title: 'Delete Notifications',
      content: `Are you sure you want to delete ${ids.length} notification(s)?`,
      okText: 'Delete',
      okType: 'danger',
      cancelText: 'Cancel',
      onOk: async () => {
        try {
          // TODO: Call API
          // await notificationService.deleteMultiple(ids)
          setNotifications(notifications.filter((n) => !ids.includes(n.id)))
          setSelectedIds([])
          message.success(`${ids.length} notification(s) deleted`)
          fetchStats()
        } catch (error) {
          message.error('Failed to delete notifications')
        }
      },
    })
  }

  const handleDeleteAllRead = async () => {
    Modal.confirm({
      title: 'Delete All Read Notifications',
      content: 'Are you sure you want to delete all read notifications?',
      okText: 'Delete',
      okType: 'danger',
      cancelText: 'Cancel',
      onOk: async () => {
        try {
          // TODO: Call API
          // await notificationService.deleteAllRead()
          setNotifications(notifications.filter((n) => !n.read))
          message.success('All read notifications deleted')
          fetchStats()
        } catch (error) {
          message.error('Failed to delete read notifications')
        }
      },
    })
  }

  const handleNotificationClick = (notification: Notification) => {
    if (!notification.read) {
      handleMarkAsRead([notification.id])
    }
    if (notification.actionUrl) {
      navigate(notification.actionUrl)
    }
  }

  const getNotificationIcon = (type: NotificationType) => {
    switch (type) {
      case 'success':
        return <CheckCircleOutlined style={{ fontSize: 24, color: DesignTokens.colors.success.main }} />
      case 'error':
        return <CloseCircleOutlined style={{ fontSize: 24, color: DesignTokens.colors.error.main }} />
      case 'warning':
        return <WarningOutlined style={{ fontSize: 24, color: DesignTokens.colors.warning.main }} />
      default:
        return <InfoCircleOutlined style={{ fontSize: 24, color: DesignTokens.colors.info.main }} />
    }
  }

  const getCategoryColor = (category: NotificationCategory): string => {
    const colors: Record<NotificationCategory, string> = {
      pipeline: 'blue',
      task: 'purple',
      system: 'cyan',
      security: 'red',
      project: 'green',
      sample: 'orange',
      result: 'magenta',
    }
    return colors[category] || 'default'
  }

  const getTimeAgo = (timestamp: string): string => {
    const date = new Date(timestamp)
    const now = new Date()
    const diffMs = now.getTime() - date.getTime()
    const diffMins = Math.floor(diffMs / 60000)
    const diffHours = Math.floor(diffMs / 3600000)
    const diffDays = Math.floor(diffMs / 86400000)

    if (diffMins < 1) return 'Just now'
    if (diffMins < 60) return `${diffMins} minutes ago`
    if (diffHours < 24) return `${diffHours} hours ago`
    if (diffDays === 1) return 'Yesterday'
    if (diffDays < 7) return `${diffDays} days ago`
    return date.toLocaleDateString()
  }

  return (
    <div className="notifications-page">
      <PageHeader
        title="通知中心"
        subtitle="查看和管理您的所有通知"
        icon={<BellOutlined />}
      />

      {/* Statistics */}
      <Row gutter={[16, 16]} style={{ marginBottom: 24 }}>
        <Col xs={12} sm={6}>
          <StatisticCard
            title="全部通知"
            value={stats.total}
            icon={<BellOutlined />}
            color={DesignTokens.colors.primary.main}
          />
        </Col>
        <Col xs={12} sm={6}>
          <StatisticCard
            title="未读通知"
            value={stats.unread}
            icon={<InboxOutlined />}
            color={DesignTokens.colors.warning.main}
          />
        </Col>
        <Col xs={12} sm={6}>
          <StatisticCard
            title="成功通知"
            value={stats.success}
            icon={<CheckCircleOutlined />}
            color={DesignTokens.colors.success.main}
          />
        </Col>
        <Col xs={12} sm={6}>
          <StatisticCard
            title="错误通知"
            value={stats.error}
            icon={<CloseCircleOutlined />}
            color={DesignTokens.colors.error.main}
          />
        </Col>
      </Row>

      {/* Filters and Actions */}
      <Card className="notifications-card">
        <div className="notifications-toolbar">
          <Space wrap>
            <Checkbox
              checked={selectedIds.length === notifications.length && notifications.length > 0}
              indeterminate={
                selectedIds.length > 0 && selectedIds.length < notifications.length
              }
              onChange={(e) => handleSelectAll(e.target.checked)}
            >
              Select All
            </Checkbox>

            {selectedIds.length > 0 && (
              <>
                <Button
                  icon={<CheckOutlined />}
                  onClick={() => handleMarkAsRead(selectedIds)}
                >
                  Mark as Read ({selectedIds.length})
                </Button>
                <Button
                  danger
                  icon={<DeleteOutlined />}
                  onClick={() => handleDelete(selectedIds)}
                >
                  Delete ({selectedIds.length})
                </Button>
              </>
            )}

            {selectedIds.length === 0 && (
              <>
                <Button icon={<CheckOutlined />} onClick={handleMarkAllAsRead}>
                  Mark All as Read
                </Button>
                <Button danger icon={<DeleteOutlined />} onClick={handleDeleteAllRead}>
                  Delete All Read
                </Button>
              </>
            )}
          </Space>

          <Space wrap>
            <Select
              value={filterType}
              onChange={setFilterType}
              style={{ width: 120 }}
              suffixIcon={<FilterOutlined />}
            >
              <Option value="all">All Types</Option>
              <Option value="success">Success</Option>
              <Option value="error">Error</Option>
              <Option value="warning">Warning</Option>
              <Option value="info">Info</Option>
            </Select>

            <Select
              value={filterCategory}
              onChange={setFilterCategory}
              style={{ width: 120 }}
              suffixIcon={<FilterOutlined />}
            >
              <Option value="all">All Categories</Option>
              <Option value="pipeline">Pipeline</Option>
              <Option value="task">Task</Option>
              <Option value="system">System</Option>
              <Option value="security">Security</Option>
              <Option value="project">Project</Option>
              <Option value="sample">Sample</Option>
              <Option value="result">Result</Option>
            </Select>

            <Select
              value={filterRead}
              onChange={setFilterRead}
              style={{ width: 100 }}
              suffixIcon={<FilterOutlined />}
            >
              <Option value="all">All</Option>
              <Option value="unread">Unread</Option>
              <Option value="read">Read</Option>
            </Select>
          </Space>
        </div>

        {/* Notifications List */}
        {loading ? (
          <div style={{ textAlign: 'center', padding: '48px' }}>
            <Spin size="large" />
          </div>
        ) : notifications.length === 0 ? (
          <Empty
            image={Empty.PRESENTED_IMAGE_SIMPLE}
            description="No notifications found"
            style={{ padding: '48px' }}
          />
        ) : (
          <>
            <List
              dataSource={notifications}
              renderItem={(item) => (
                <List.Item
                  className={`notification-list-item ${!item.read ? 'unread' : ''}`}
                  onClick={() => handleNotificationClick(item)}
                  actions={[
                    !item.read && (
                      <Button
                        type="text"
                        size="small"
                        icon={<CheckOutlined />}
                        onClick={(e) => {
                          e.stopPropagation()
                          handleMarkAsRead([item.id])
                        }}
                      >
                        Mark as Read
                      </Button>
                    ),
                    item.actionUrl && (
                      <Button
                        type="link"
                        size="small"
                        onClick={(e) => {
                          e.stopPropagation()
                          handleNotificationClick(item)
                        }}
                      >
                        {item.actionText || 'View'}
                      </Button>
                    ),
                    <Button
                      type="text"
                      size="small"
                      danger
                      icon={<DeleteOutlined />}
                      onClick={(e) => {
                        e.stopPropagation()
                        handleDelete([item.id])
                      }}
                    />,
                  ].filter(Boolean)}
                >
                  <div className="notification-select">
                    <Checkbox
                      checked={selectedIds.includes(item.id)}
                      onChange={(e) => {
                        e.stopPropagation()
                        handleSelect(item.id, e.target.checked)
                      }}
                      onClick={(e) => e.stopPropagation()}
                    />
                  </div>

                  <List.Item.Meta
                    avatar={getNotificationIcon(item.type)}
                    title={
                      <div className="notification-list-title">
                        <Space>
                          {item.title}
                          {!item.read && <span className="unread-badge">NEW</span>}
                        </Space>
                        <Tag color={getCategoryColor(item.category)}>
                          {item.category.toUpperCase()}
                        </Tag>
                      </div>
                    }
                    description={
                      <div className="notification-list-description">
                        <Paragraph type="secondary" style={{ marginBottom: 4 }}>
                          {item.message}
                        </Paragraph>
                        <Text type="secondary" style={{ fontSize: 12 }}>
                          {getTimeAgo(item.timestamp)}
                        </Text>
                      </div>
                    }
                  />
                </List.Item>
              )}
            />

            {/* Pagination */}
            <div className="notifications-pagination">
              <Pagination
                current={currentPage}
                pageSize={pageSize}
                total={total}
                onChange={(page, size) => {
                  setCurrentPage(page)
                  setPageSize(size || 10)
                }}
                showSizeChanger
                showTotal={(total) => `Total ${total} notifications`}
              />
            </div>
          </>
        )}
      </Card>
    </div>
  )
}
